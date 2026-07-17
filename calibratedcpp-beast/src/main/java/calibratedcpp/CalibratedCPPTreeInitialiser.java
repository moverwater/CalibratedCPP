package calibratedcpp;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.util.Randomizer;
import calibration.CalibrationForest;
import calibration.CalibrationNode;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import org.hipparchus.distribution.continuous.BetaDistribution;
import org.hipparchus.special.Erf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Marcus Overwater
 */
@Description("Simulates a starting tree from a calibrated coalescent point process whose "
        + "calibration clades are monophyletic and within the support of their MRCA priors.")
public class CalibratedCPPTreeInitialiser extends Tree implements StateNodeInitialiser {

    public Input<CalibratedCoalescentPointProcess> modelInput =
            new Input<>("model", "The calibrated CPP tree distribution supplying the node-age law "
                    + "(via its CDF), the origin/root conditioning, and — when no clade priors are "
                    + "given — the calibration clade structure. Optional: when omitted the initialiser "
                    + "discovers the calibrated-CPP prior (and the MRCA/Calibration priors) attached "
                    + "to the tree, so the model stays defined only in the posterior.");

    public Input<List<MRCAPrior>> constraintInput =
            new Input<>("constraint", "MRCA priors constraining clade ages. Each clade's starting "
                    + "MRCA age is drawn within the support of its distribution (e.g. Uniform "
                    + "bounds, or an offset/truncation), honouring monophyly.", new ArrayList<>());

    public Input<CalibrationForest> calibrationForestInput =
            new Input<>("calibrationForest", "Calibration clades with soft age bounds "
                    + "(the CalibrationPrior path). Alternative/complement to 'constraint'.");

    public Input<Distribution> calibrationInput =
            new Input<>("calibration", "A calibration prior (or a CompoundDistribution such as "
                    + "CalibrationDistribution wrapping one) to read clade bounds from. MRCAPrior "
                    + "children contribute hard support; a CalibrationPrior contributes its soft "
                    + "clade bounds. Lets BEAUti wire the calibration panel's output with one idref.");

    public Input<Double> rootHeightInput =
            new Input<>("rootHeight", "Optional starting root age. Used when conditioning on the "
                    + "root; otherwise the root age is sampled below the origin.");

    // Smallest gap enforced between a parent age and its oldest child, so heights are strictly
    // increasing towards the root even after clamping to feasible bounds.
    private static final double MIN_GAP = 1e-8;

    private CalibratedCoalescentPointProcess model;
    private List<String> taxaNames;
    private Set<String> realTaxaSet;

    // Per-clade MRCA-age prior, keyed by the clade's (order-independent) taxon-id set.
    private Map<Set<String>, CladeAge> cladeBounds;

    /**
     * How to draw a clade's starting MRCA age: its support [lower, upper] plus the prior it is drawn
     * from — an {@link MRCAPrior}'s distribution (absolute age, truncated to the nesting window), or a
     * {@link CalibrationPrior}'s fitted hierarchy (top clade: LogNormal(μ,σ) absolute age; nested +
     * overlapping: Beta(α,β) on the child/parent age ratio; nested + non-overlapping: LogNormal
     * truncated at the parent age). With neither, the age is drawn from the CPP node-age law.
     */
    private static final class CladeAge {
        double lower = 0.0;
        double upper = Double.POSITIVE_INFINITY;
        ScalarDistribution<?, ?> mrcaDistr;   // MRCAPrior path (null otherwise)
        boolean calibration;                  // CalibrationPrior path
        double mu, sigma;                     // fitted LogNormal (log-mean, log-sd)
        boolean overlapEdge;                  // nested + overlapping parent → Beta on the ratio
        double alpha, beta;                   // fitted Beta parameters

        static CladeAge mrca(double lower, double upper, ScalarDistribution<?, ?> distr) {
            CladeAge a = new CladeAge();
            a.lower = lower; a.upper = upper; a.mrcaDistr = distr;
            return a;
        }

        static CladeAge calibration(CalibrationCladePrior p) {
            CladeAge a = new CladeAge();
            a.lower = p.getLower(); a.upper = p.getUpper();
            if (p.sigma2 > 0.0) { // fitted by CalibrationPrior.initAndValidate()
                a.calibration = true;
                a.mu = p.mu; a.sigma = Math.sqrt(p.sigma2);
                a.overlapEdge = p.isOverlapEdge; a.alpha = p.alpha; a.beta = p.beta;
            }
            return a;
        }
    }

    @Override
    public void initAndValidate() {
        taxaNames = new ArrayList<>();
        if (m_taxonset.get() != null) {
            taxaNames.addAll(m_taxonset.get().asStringList());
        } else if (m_initial.get() != null) {
            for (Node leaf : m_initial.get().getExternalNodes()) {
                taxaNames.add(leaf.getID());
            }
        } else {
            throw new IllegalArgumentException(
                    "CalibratedCPPTreeInitialiser needs a 'taxonset' or an 'initial' tree to "
                    + "know the taxa.");
        }
        if (taxaNames.isEmpty()) {
            throw new IllegalArgumentException("CalibratedCPPTreeInitialiser: no taxa found.");
        }
        realTaxaSet = new LinkedHashSet<>(taxaNames);

        initStateNodes();
        super.initAndValidate();
    }

    @Override
    public void initStateNodes() {

        model = resolveModel();
        if (model == null) {
            return;
        }
        // Readiness gate: updateModel() builds the node-age law and throws if the model's rate
        // parameters are not wired yet (as during BEAUti's createSubNet). In that case do nothing.
        try {
            model.updateModel();
        } catch (Throwable e) {
            return;
        }

        Node root = simulate();

        final int nTaxa = taxaNames.size();
        leafNodeCount = nTaxa;
        internalNodeCount = nTaxa - 1;
        nodeCount = 2 * nTaxa - 1;

        // Preserve the initial tree's taxon numbering where possible, exactly as RandomTree does,
        // so the tree we write in lines up with operators/loggers referencing node numbers.
        Map<String, Integer> taxonToNr = new HashMap<>();
        if (m_initial.get() != null && m_initial.get().getLeafNodeCount() == nTaxa) {
            for (Node n : m_initial.get().getExternalNodes()) {
                taxonToNr.put(n.getID(), n.getNr());
            }
        } else {
            for (int k = 0; k < taxaNames.size(); k++) {
                taxonToNr.put(taxaNames.get(k), k);
            }
        }

        numberNodes(root, new int[]{0}, taxonToNr);

        // Assign the root field directly and let initArrays() (re)allocate m_nodes to nodeCount, as
        // RandomTree does. setRoot() must NOT be used here: it indexes the *existing* m_nodes by the
        // new root's number, and when this initialiser no-opped at parse (model not yet discovered)
        // super.initAndValidate left m_nodes as a length-1 dummy → ArrayIndexOutOfBounds at run time.
        this.root = root;
        initArrays();

        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        if (m_initial.get() != null) {
            stateNodes.add(m_initial.get());
        }
    }

    /**
     * The calibrated-CPP tree prior to initialise from: the explicit {@code model} input if given,
     * otherwise the (first) {@link CalibratedCoalescentPointProcess} attached to the target tree.
     * Discovering it from the tree lets the generated XML leave the model defined only in the
     * posterior, with the {@code <init>} element carrying no model reference at all.
     */
    private CalibratedCoalescentPointProcess resolveModel() {
        if (modelInput.get() != null) {
            return modelInput.get();
        }
        Tree tree = m_initial.get();
        if (tree == null) {
            return null;
        }
        for (BEASTInterface out : tree.getOutputs()) {
            if (out instanceof CalibratedCoalescentPointProcess cpp) {
                return cpp;
            }
        }
        return null;
    }

    // ------------------------------------------------------------------
    // Simulation

    private Node simulate() {
        // Gather clade structure (taxon sets) and their age bounds from the priors in play. Clade
        // taxa are resolved to the tree's real taxa first (see collectCladeBounds / resolveId).
        // (The node-age law was already built by the updateModel() readiness gate in initStateNodes.)
        List<TaxonSet> cladeTaxonSets = collectCladeBounds();

        Node root = null;
        if (!cladeTaxonSets.isEmpty()) {
            try {
                root = buildFromClades(cladeTaxonSets);
            } catch (RuntimeException e) {
                Log.warning.println("CalibratedCPPTreeInitialiser: could not build a clade-constrained "
                        + "starting tree (" + e.getMessage() + "); falling back to an unconstrained tree.");
                root = null;
            }
        }

        // Safety net: the written tree must contain exactly the tree's taxa, each exactly once
        // (checked by count AND set, so a duplicated-plus-missing leaf can't slip through and later
        // collide node numbers in numberNodes/initArrays). Otherwise fall back to a plain coalescent.
        if (root == null || !coversAllTaxaExactly(root)) {
            if (root != null) {
                Log.warning.println("CalibratedCPPTreeInitialiser: clade-based tree did not cover all "
                        + "taxa exactly once; falling back to an unconstrained starting tree.");
            }
            root = plainTree();
        }
        enforceHeights(root);
        return root;
    }

    /** Build a starting tree honouring the calibration clades (already resolved to real taxa). */
    private Node buildFromClades(List<TaxonSet> cladeTaxonSets) {
        final int nTaxa = taxaNames.size();
        CalibrationForest forest = new CalibrationForest(cladeTaxonSets);

        // Detect a whole-tree (root) clade; the model strips it, but a prior may still target it.
        CalibrationNode rootClade = null;
        for (CalibrationNode node : forest.getRoots()) {
            if (taxaOf(node).size() == nTaxa) {
                rootClade = node;
                break;
            }
        }

        List<CalibrationNode> topClades =
                (rootClade != null) ? rootClade.children : forest.getRoots();

        // Floor for the root: it must sit above every calibration clade's required floor.
        double topFloor = 0.0;
        for (CalibrationNode c : topClades) topFloor = Math.max(topFloor, requiredFloor(c));

        // conditionAge is the origin (stem) when conditioning on it, else the root age. The tree is
        // built conditioned on that age, exactly as the LPhy sampler conditions on stem or root.
        boolean rootConditioned = model.isConditionOnRoot() || model.getOrigin() == null;
        double conditionAge;
        Map<CalibrationNode, Double> ages = new HashMap<>();
        if (!rootConditioned) {
            conditionAge = model.getOrigin();
        } else if (rootClade != null) {
            // The whole-tree clade IS the root: draw its age from its own prior (uncapped above, so it
            // isn't clamped below its calibration by an arbitrary chooseRootAge). Only fall back to
            // chooseRootAge when the root clade carries no prior at all.
            double rootCap = cladeBounds.containsKey(taxaOf(rootClade))
                    ? Double.POSITIVE_INFINITY : chooseRootAge(topFloor);
            conditionAge = sampleCladeAge(topFloor, rootCap, rootClade);
            ages.put(rootClade, conditionAge);
        } else {
            conditionAge = chooseRootAge(topFloor);
        }

        // Choose each clade's MRCA age within its prior support, respecting nesting (child < parent).
        for (CalibrationNode c : topClades) assignAges(c, conditionAge, ages);

        return buildCPPSubtree(new ArrayList<>(taxaNames), topClades, conditionAge, rootConditioned, ages);
    }

    /** Unconstrained fallback: a plain CPP over all taxa (the no-clade case of {@link #buildCPPSubtree}). */
    private Node plainTree() {
        boolean rootConditioned = model.isConditionOnRoot() || model.getOrigin() == null;
        double conditionAge = rootConditioned ? chooseRootAge(0.0) : model.getOrigin();
        return buildCPPSubtree(new ArrayList<>(taxaNames), Collections.emptyList(),
                conditionAge, rootConditioned, new HashMap<>());
    }

    /**
     * Port of {@code AbstractCalibratedCPPTree.sample()}: build a CPP (sub)tree over {@code taxa}
     * conditioned on {@code conditionAge}, with {@code clades} the maximal calibrations at this level
     * (each carrying a chosen age in {@code ages}). Maximal clades are dropped into weighted "slots"
     * — the weight of a slot is {@code (Q(conditionAge) - Q(cladeAge))^s} where {@code s} counts its
     * still-free neighbouring node-age positions — their subtrees are built recursively, remaining
     * node ages are drawn from the law truncated below {@code conditionAge}, and the linear
     * arrangement is coalesced youngest-first. When conditioning on the root, {@code conditionAge}
     * is planted at a random slot so it becomes the root height; otherwise it is the (stem) origin
     * and the root falls below it.
     */
    private Node buildCPPSubtree(List<String> taxa, List<CalibrationNode> clades,
                                 double conditionAge, boolean rootConditioned,
                                 Map<CalibrationNode, Double> ages) {
        // Process maximal clades oldest-first (LPhy sorts the calibrations by descending age).
        List<CalibrationNode> sorted = new ArrayList<>(clades);
        sorted.sort((a, b) -> Double.compare(ages.get(b), ages.get(a)));

        Set<String> cladeTaxa = new HashSet<>();
        for (CalibrationNode c : clades) cladeTaxa.addAll(taxaOf(c));
        List<String> freeTaxa = new ArrayList<>();
        for (String t : taxa) if (!cladeTaxa.contains(t)) freeTaxa.add(t);

        // m = free tips + one slot per maximal clade (each clade collapses to its subtree root).
        int m = freeTaxa.size() + clades.size();
        if (m <= 1) {
            return clades.isEmpty() ? leaf(freeTaxa.get(0)) : buildClade(clades.get(0), ages);
        }

        List<Integer> free = new ArrayList<>();       // available slot indices (LPhy's "A")
        for (int i = 0; i < m; i++) free.add(i);
        List<Double> times = new ArrayList<>(Collections.nCopies(m, 0.0));
        List<Node> nodeList = new ArrayList<>(Collections.nCopies(m, (Node) null));

        // Plant the root age at a random interior slot so it ends up as the root height.
        if (rootConditioned) {
            int ind = (m == 2) ? 1 : Randomizer.nextInt(m - 1) + 1;
            times.set(ind, conditionAge);
        }

        for (CalibrationNode clade : sorted) {
            double cladeAge = ages.get(clade);
            double qCond = cdf(conditionAge);
            double qClade = cdf(cladeAge);
            double w = (Double.isNaN(qCond) || Double.isNaN(qClade)) ? 1.0 : (qCond - qClade);
            if (!(w > 0.0)) w = 1e-12; // keep slot weights positive/finite even at degenerate ages

            int[] s = calculateScore(free, m, times);
            double[] weights = getWeights(s, w);
            int li = (free.size() == 1) ? free.get(0) : free.get(sampleIndex(weights));

            nodeList.set(li, buildClade(clade, ages));
            if (li > 0 && times.get(li) == 0) times.set(li, sampleAge(cladeAge, conditionAge));
            if (li < m - 1 && times.get(li + 1) == 0) times.set(li + 1, sampleAge(cladeAge, conditionAge));
            free.remove(Integer.valueOf(li));
        }

        // Fill the remaining node-age positions from the law truncated to (0, conditionAge).
        for (int i = 1; i < m; i++) {
            if (times.get(i) == 0) times.set(i, sampleAge(0.0, conditionAge));
        }

        // Slot 0 is never a merge point; set it to the maximum so it can't be picked as "youngest".
        if (rootConditioned) {
            double max = times.get(1);
            for (int i = 2; i < m; i++) max = Math.max(max, times.get(i));
            times.set(0, max);
        } else {
            times.set(0, conditionAge);
        }

        // Drop the uncalibrated taxa into the empty slots in random order.
        shuffle(freeTaxa);
        int idx = 0;
        for (int i = 0; i < m && idx < freeTaxa.size(); i++) {
            if (nodeList.get(i) == null) nodeList.set(i, leaf(freeTaxa.get(idx++)));
        }

        coalesce(nodeList, times);
        return nodeList.get(0);
    }

    /**
     * Collects the calibration clades and their MRCA-age bounds from {@code constraint} MRCA priors
     * and any {@code calibrationForest}, falling back to the model's (bound-free) calibration clades.
     * Populates {@link #cladeBounds} and returns the deduplicated clade taxon sets.
     */
    private List<TaxonSet> collectCladeBounds() {
        cladeBounds = new HashMap<>();
        List<TaxonSet> taxonSets = new ArrayList<>();
        Set<Set<String>> seen = new HashSet<>();

        // Soft bounds from a CalibrationForest (CalibrationPrior path).
        CalibrationForest forest = calibrationForestInput.get();
        if (forest != null) {
            for (CalibrationNode node : forest.getAllNodes()) {
                addClade(node.taxa, node.prior != null ? CladeAge.calibration(node.prior) : null,
                        taxonSets, seen);
            }
        }

        // Bounds unpacked from a calibration prior / CompoundDistribution wrapper (the BEAUti panel).
        if (calibrationInput.get() != null) {
            unpackCalibration(calibrationInput.get(), taxonSets, seen);
        }

        // Hard support from MRCA priors — takes precedence, as it reflects what the sampler enforces.
        for (MRCAPrior prior : constraintInput.get()) {
            addMRCAPrior(prior, taxonSets, seen);
        }

        // Runtime discovery: when nothing was passed explicitly, read the MRCA/Calibration priors
        // attached to the tree. This lets the generated XML leave the calibration priors defined in
        // the posterior, with the <init> element referencing none of them.
        if (taxonSets.isEmpty()) {
            Tree tree = m_initial.get();
            if (tree != null) {
                for (BEASTInterface out : tree.getOutputs()) {
                    if (out instanceof MRCAPrior mp) {
                        addMRCAPrior(mp, taxonSets, seen);
                    } else if (out instanceof CalibrationPrior cp) {
                        unpackCalibration(cp, taxonSets, seen);
                    }
                }
            }
        }

        // If nothing carried bounds, at least take clade structure from the model's calibrations.
        if (taxonSets.isEmpty() && model.getCalibrations() != null) {
            for (TaxonSet ts : model.getCalibrations()) {
                addClade(ts, null, taxonSets, seen);
            }
        }

        return taxonSets;
    }

    /**
     * Recursively reads clade bounds from a calibration prior: a {@link CompoundDistribution}
     * (e.g. the {@code CalibrationDistribution} wrapper) is descended into; an {@link MRCAPrior}
     * contributes the hard support of its distribution; a {@link CalibrationPrior} contributes its
     * soft clade bounds (from its explicit clades or its calibration forest).
     */
    private void unpackCalibration(Distribution dist, List<TaxonSet> taxonSets, Set<Set<String>> seen) {
        if (dist instanceof CompoundDistribution compound) {
            for (Distribution child : compound.pDistributions.get()) {
                unpackCalibration(child, taxonSets, seen);
            }
        } else if (dist instanceof MRCAPrior prior) {
            addMRCAPrior(prior, taxonSets, seen);
        } else if (dist instanceof CalibrationPrior cp) {
            for (CalibrationCladePrior clade : cp.cladesInput.get()) {
                addClade(clade.getTaxa(), CladeAge.calibration(clade), taxonSets, seen);
            }
            CalibrationForest cpForest = cp.calibrationForestInput.get();
            if (cpForest != null) {
                for (CalibrationNode node : cpForest.getAllNodes()) {
                    addClade(node.taxa, node.prior != null ? CladeAge.calibration(node.prior) : null,
                            taxonSets, seen);
                }
            }
        }
    }

    private void addMRCAPrior(MRCAPrior prior, List<TaxonSet> taxonSets, Set<Set<String>> seen) {
        TaxonSet ts = prior.taxonsetInput.get();
        if (ts == null) return;
        double[] sb = supportBounds(prior);
        addClade(ts, CladeAge.mrca(sb[0], sb[1], prior.distInput.get()), taxonSets, seen);
    }

    /**
     * Register a clade (once) and its age prior (if any), keyed by its resolved real-taxon-id set.
     * Clade taxa are first mapped to the tree's real taxa via {@link #resolveId(String)} — this
     * repairs BEAUti's duplicate suffixed taxa (e.g. {@code Homo_sapiens1} → {@code Homo_sapiens})
     * so nesting is recovered and no phantom leaves ever enter the tree. Clades that resolve to
     * fewer than two real taxa carry no internal-node age and are skipped.
     */
    private void addClade(TaxonSet rawTaxa, CladeAge age, List<TaxonSet> taxonSets, Set<Set<String>> seen) {
        Set<String> real = resolveTaxa(taxaIds(rawTaxa));
        if (real.size() < 2) return;
        if (seen.add(real)) taxonSets.add(realTaxonSet(real));
        if (age != null) cladeBounds.put(real, age);
    }

    /** Map a clade taxon id to a real tree taxon, repairing BEAUti's duplicate suffixes; null if none. */
    private String resolveId(String id) {
        if (realTaxaSet.contains(id)) return id;
        String base = id.replaceAll("[0-9]+$", "");
        return (!base.equals(id) && realTaxaSet.contains(base)) ? base : null;
    }

    private Set<String> resolveTaxa(Set<String> ids) {
        Set<String> out = new LinkedHashSet<>();
        for (String id : ids) {
            String r = resolveId(id);
            if (r != null) out.add(r);
        }
        return out;
    }

    /** A TaxonSet of real taxa (id-based), used as clade structure for the forest builder. */
    private TaxonSet realTaxonSet(Set<String> realIds) {
        List<Taxon> taxa = new ArrayList<>();
        for (String id : realIds) taxa.add(new Taxon(id));
        return new TaxonSet(taxa);
    }

    /** True iff {@code root}'s leaves are exactly the tree's real taxa, each present exactly once. */
    private boolean coversAllTaxaExactly(Node root) {
        List<String> ids = new ArrayList<>();
        collectLeafIds(root, ids);
        return ids.size() == realTaxaSet.size() && new HashSet<>(ids).equals(realTaxaSet);
    }

    private void collectLeafIds(Node node, List<String> ids) {
        if (node.isLeaf()) {
            ids.add(node.getID());
        } else {
            for (Node child : node.getChildren()) collectLeafIds(child, ids);
        }
    }

    /** Hard [lower, upper] support of an MRCA prior's distribution; [0, +inf] if monophyly-only. */
    private double[] supportBounds(MRCAPrior prior) {
        ScalarDistribution<?, ?> distr = prior.distInput.get();
        if (distr == null) {
            return new double[]{0.0, Double.POSITIVE_INFINITY};
        }
        double lower = 0.0;
        double upper = Double.POSITIVE_INFINITY;
        try {
            lower = ((Number) distr.inverseCumulativeProbability(0.0)).doubleValue();
            upper = ((Number) distr.inverseCumulativeProbability(1.0)).doubleValue();
            if (lower > upper) {
                double tmp = lower;
                lower = upper;
                upper = tmp;
            }
        } catch (Exception e) {
            Log.warning.println("CalibratedCPPTreeInitialiser: could not read bounds from MRCAPrior "
                    + prior.getID() + " (" + e.getMessage() + "); treating as unbounded.");
        }
        if (Double.isNaN(lower) || Double.isInfinite(lower)) lower = 0.0;
        if (Double.isNaN(upper)) upper = Double.POSITIVE_INFINITY;
        return new double[]{Math.max(0.0, lower), upper};
    }

    /**
     * Build a clade's subtree as a root-conditioned CPP over the clade's taxa, with the clade's
     * direct children as the nested maximal calibrations and the clade's chosen age as the root —
     * the recursive {@code newSubClade} step of the LPhy sampler.
     */
    private Node buildClade(CalibrationNode clade, Map<CalibrationNode, Double> ages) {
        List<String> taxa = new ArrayList<>(taxaOf(clade));
        if (taxa.size() == 1) {
            return leaf(taxa.get(0)); // degenerate single-taxon clade: no MRCA node
        }
        return buildCPPSubtree(taxa, clade.children, ages.get(clade), true, ages);
    }

    /** Assign an age to {@code clade} strictly below {@code cap}, then recurse into its children. */
    private void assignAges(CalibrationNode clade, double cap, Map<CalibrationNode, Double> ages) {
        double lo = requiredFloor(clade);
        double age = sampleCladeAge(lo, cap, clade);
        ages.put(clade, age);
        for (CalibrationNode child : clade.children) {
            assignAges(child, age, ages);
        }
    }

    /**
     * Highest lower bound that {@code clade}'s MRCA must exceed for its whole subtree to remain
     * feasible: the max of its own lower bound and the required floors of all descendant clades.
     */
    private double requiredFloor(CalibrationNode clade) {
        double floor = bounds(clade)[0];
        for (CalibrationNode child : clade.children) {
            floor = Math.max(floor, requiredFloor(child));
        }
        return floor;
    }

    /**
     * Draw {@code clade}'s starting MRCA age from its prior, constrained to the nesting window
     * {@code (floor, cap)} where {@code cap} is the parent's age. MRCAPriors draw from their own
     * distribution truncated to the window; a CalibrationPrior draws hierarchically (top clade:
     * LogNormal absolute age; nested overlapping: Beta on the child/parent ratio × the parent age;
     * nested non-overlapping: LogNormal truncated at the parent). With no usable prior, the age comes
     * from the CPP node-age law. The result is clamped into the window (and a warning is issued when
     * the window is empty, i.e. the calibration bounds conflict with the nesting).
     */
    private double sampleCladeAge(double floor, double cap, CalibrationNode clade) {
        CladeAge a = cladeBounds.get(taxaOf(clade));
        double lo = Math.max(floor, (a != null) ? a.lower : 0.0);

        double t = Double.NaN;
        if (a != null && a.mrcaDistr != null) {
            t = truncatedDistrSample(a.mrcaDistr, lo, cap);
        } else if (a != null && a.calibration) {
            boolean nested = clade.parent != null && !Double.isInfinite(cap);
            if (nested && a.overlapEdge && a.alpha > 0.0 && a.beta > 0.0) {
                t = betaSample(a.alpha, a.beta) * cap;              // relative age × parent age
            } else if (nested) {
                t = truncatedLogNormalSample(a.mu, a.sigma, cap);  // truncated at the parent age
            } else {
                t = logNormalSample(a.mu, a.sigma);                // top-level: absolute age
            }
        } else if (a != null && !Double.isInfinite(a.upper)) {
            // A bounded interval with no fitted distribution (e.g. a CalibrationCladePrior not run
            // through a CalibrationPrior): draw uniformly within [lo, min(cap, upper)].
            double hi = Math.min(cap, a.upper);
            if (hi > lo) t = lo + Randomizer.nextDouble() * (hi - lo);
        }
        if (Double.isNaN(t) || Double.isInfinite(t)) {
            t = (cap > lo) ? sampleAge(lo, cap) : cap;             // node-age law / degenerate fallback
        }
        return clampToWindow(t, lo, cap, clade, a);
    }

    private double clampToWindow(double t, double lo, double cap, CalibrationNode clade, CladeAge a) {
        if (!(cap > lo)) {
            if (a != null) {
                Log.warning.println("CalibratedCPPTreeInitialiser: calibration clade " + cladeLabel(clade)
                        + " has no feasible MRCA age — its bounds [" + a.lower + ", " + a.upper + "] "
                        + "conflict with the nesting constraints ([" + lo + ", " + cap + "]). The starting "
                        + "age is clamped and will violate the prior; please check the calibration bounds.");
            }
            return Math.max(MIN_GAP, cap);
        }
        return Math.min(Math.max(t, lo), cap);
    }

    /** LogNormal(μ, σ) draw. */
    private static double logNormalSample(double mu, double sigma) {
        return Math.exp(mu + sigma * Randomizer.nextGaussian());
    }

    /** LogNormal(μ, σ) draw truncated to (0, upper) by inverse-CDF (matches CalibrationPrior's Erf math). */
    private static double truncatedLogNormalSample(double mu, double sigma, double upper) {
        double cdfUpper = 0.5 * (1.0 + Erf.erf((Math.log(upper) - mu) / (sigma * Math.sqrt(2.0))));
        if (!(cdfUpper > 0.0)) return upper * 0.5; // parent below the whole distribution: degenerate
        double p = Randomizer.nextDouble() * cdfUpper;
        return Math.exp(mu + sigma * Math.sqrt(2.0) * Erf.erfInv(2.0 * p - 1.0));
    }

    /** Beta(α, β) draw by inverse-CDF; NaN if the solver fails (caller falls back to the node-age law). */
    private static double betaSample(double alpha, double beta) {
        try {
            return new BetaDistribution(alpha, beta).inverseCumulativeProbability(Randomizer.nextDouble());
        } catch (Exception e) {
            return Double.NaN;
        }
    }

    /** Draw from a scalar distribution truncated to (lo, cap); NaN if the window is empty/unusable. */
    private static double truncatedDistrSample(ScalarDistribution<?, ?> distr, double lo, double cap) {
        try {
            double cLo = (lo <= 0.0) ? 0.0 : distr.cumulativeProbability(lo);
            double cHi = Double.isInfinite(cap) ? 1.0 : distr.cumulativeProbability(cap);
            if (!(cHi > cLo)) return Double.NaN;
            double u = cLo + Randomizer.nextDouble() * (cHi - cLo);
            return ((Number) distr.inverseCumulativeProbability(u)).doubleValue();
        } catch (Exception e) {
            return Double.NaN;
        }
    }

    /** [lower, upper] MRCA-age bounds for a clade, defaulting to [0, +inf] when unconstrained. */
    private double[] bounds(CalibrationNode clade) {
        CladeAge a = cladeBounds.get(taxaOf(clade));
        return (a != null) ? new double[]{a.lower, a.upper} : new double[]{0.0, Double.POSITIVE_INFINITY};
    }

    /** Root age: rootHeight input if given, else sampled below the origin (or above the clade floor). */
    private double chooseRootAge(double floor) {
        if (rootHeightInput.get() != null) {
            return Math.max(rootHeightInput.get(), floor + MIN_GAP);
        }
        Double origin = model.getOrigin();
        if (!model.isConditionOnRoot() && origin != null) {
            double lo = Math.min(floor, origin - MIN_GAP);
            return sampleAge(lo, origin);
        }
        // Conditioning on the root with no explicit origin: draw a root age above the clade floor from
        // the node-age law, bounded by a finite, law-informed cap so it can never run away.
        double sampled = sampleAge(floor, lawUpperCap(floor));
        return Math.max(sampled, floor + MIN_GAP);
    }

    /**
     * A finite upper bound for root-age sampling: the age at (about) the 0.99 quantile of the node-age
     * law above {@code floor}. Doubling stops at that quantile, at the law's saturation plateau (its
     * supremum lies below 0.99), or — if the law is unavailable — at a modest default, so we never
     * expand the bracket without bound.
     */
    private double lawUpperCap(double floor) {
        double target = Math.log(0.99);
        double hi = Math.max(1.0, floor * 2.0);
        double prev = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < 60; i++) {
            double lc = logCdf(hi);
            if (Double.isNaN(lc)) return Math.max(floor * 3.0 + 1.0, 1.0); // law unavailable
            if (lc >= target) return hi;
            if (lc - prev < 1e-6) return hi; // saturation plateau below the target quantile
            prev = lc;
            hi *= 2.0;
        }
        return hi;
    }

    /**
     * Linear coalescing of the CPP arrangement (port of {@code AbstractCalibratedCPPTree.coalesce}):
     * repeatedly merge the pair of neighbours around the youngest node age, so the largest age (slot
     * 0, never itself a merge point) becomes the root. {@code nodeList} and {@code times} are indexed
     * together and shrink in place until a single root node remains.
     */
    private void coalesce(List<Node> nodeList, List<Double> times) {
        while (nodeList.size() > 1) {
            int j = indexOfMin(times);
            // Slot 0 is the root position (the maximum age) and is never a merge point. It only
            // becomes indexOfMin's answer when all remaining times are tied (possible once clade
            // ages are drawn from their priors and clamped); force j>=1 so get(j-1) stays in range.
            if (times.size() == 2 || j == 0) j = 1;
            Node left = nodeList.get(j - 1);
            Node right = nodeList.get(j);
            Node parent = new Node();
            parent.setHeight(times.get(j));
            parent.addChild(left);
            parent.addChild(right);
            nodeList.set(j - 1, parent);
            nodeList.remove(j);
            times.remove(j);
        }
    }

    private static int indexOfMin(List<Double> times) {
        int minIndex = 0;
        double minValue = times.get(0);
        for (int i = 1; i < times.size(); i++) {
            if (times.get(i) < minValue) {
                minValue = times.get(i);
                minIndex = i;
            }
        }
        return minIndex;
    }

    /** CDF Q(t) of the node-age law, or NaN when the law is unavailable (see {@link #logCdf}). */
    private double cdf(double t) {
        double v = logCdf(t);
        return Double.isNaN(v) ? Double.NaN : Math.exp(v);
    }

    /** For each free slot, how many of its two neighbouring node-age positions are still unassigned. */
    private static int[] calculateScore(List<Integer> free, int m, List<Double> times) {
        int[] s = new int[free.size()];
        for (int j = 0; j < free.size(); j++) {
            int i = free.get(j);
            int count = 0;
            if (i < m - 1 && times.get(i + 1) == 0) count++;
            if (i > 0 && times.get(i) == 0) count++;
            s[j] = count;
        }
        return s;
    }

    /** Slot weights {@code w^s}, normalised; uniform if that is degenerate (w=0, or non-finite). */
    private static double[] getWeights(int[] s, double w) {
        double[] weights = new double[s.length];
        double sum = 0.0;
        for (int i = 0; i < s.length; i++) {
            weights[i] = Math.pow(w, s[i]);
            sum += weights[i];
        }
        if (!(sum > 0.0) || Double.isNaN(sum) || Double.isInfinite(sum)) {
            java.util.Arrays.fill(weights, 1.0 / s.length);
            return weights;
        }
        for (int i = 0; i < s.length; i++) weights[i] /= sum;
        return weights;
    }

    /** Sample an index in proportion to (already-normalised) {@code weights}. */
    private static int sampleIndex(double[] weights) {
        double u = Randomizer.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < weights.length; i++) {
            cumulative += weights[i];
            if (u <= cumulative) return i;
        }
        return weights.length - 1;
    }

    /** In-place Fisher–Yates shuffle using BEAST's seeded RNG (so runs are reproducible). */
    private static void shuffle(List<String> list) {
        for (int i = list.size() - 1; i > 0; i--) {
            Collections.swap(list, i, Randomizer.nextInt(i + 1));
        }
    }

    /** Human-readable clade label for warnings: the TaxonSet id, else its taxa. */
    private String cladeLabel(CalibrationNode clade) {
        String id = (clade.taxa != null) ? clade.taxa.getID() : null;
        return (id != null && !id.isEmpty()) ? id : taxaOf(clade).toString();
    }

    /**
     * Sample a node age in {@code (lo, hi)} from the model's node-age law by numerically inverting
     * its log-CDF (monotone increasing in age). Falls back to a uniform draw when the law is
     * degenerate over the window or not yet numerically available.
     */
    private double sampleAge(double lo, double hi) {
        if (!(hi > lo)) return Math.max(lo, hi);
        boolean finiteHi = !Double.isInfinite(hi);
        double uniformHi = finiteHi ? hi : lo + Math.max(1.0, Math.abs(lo));

        // Probe the law at both ends. logCdf returns NaN when the law is unavailable (e.g. the
        // age-dependent model's spline is not built until the first likelihood evaluation, which
        // runs after initialisation) — in which case fall back to a uniform draw.
        double logQlo = logCdf(lo);
        double logQhi = logCdf(uniformHi);
        double qLo = (lo <= 0.0) ? 0.0 : Math.exp(logQlo);
        double qHi = finiteHi ? Math.exp(logQhi) : 1.0;

        if (Double.isNaN(qLo) || Double.isNaN(qHi) || !(qHi > qLo)) {
            return lo + Randomizer.nextDouble() * (uniformHi - lo);
        }

        double target = Math.log(qLo + Randomizer.nextDouble() * (qHi - qLo));
        double a = lo;
        double b = finiteHi ? hi : uniformHi;
        // Expand the upper bracket if hi is infinite until the CDF exceeds the target.
        if (!finiteHi) {
            for (int i = 0; i < 100 && logCdf(b) < target; i++) {
                a = b;
                b *= 2.0;
            }
        }
        double tol = 1e-9 * Math.max(1.0, b);
        for (int it = 0; it < 100 && (b - a) > tol; it++) {
            double mid = 0.5 * (a + b);
            if (logCdf(mid) < target) {
                a = mid;
            } else {
                b = mid;
            }
        }
        return 0.5 * (a + b);
    }

    /** Model log-CDF of the node-age law, or {@code NaN} if the law cannot be evaluated yet. */
    private double logCdf(double t) {
        try {
            double v = model.calculateLogNodeAgeCDF(t);
            return Double.isNaN(v) ? Double.NaN : v;
        } catch (Throwable e) {
            return Double.NaN;
        }
    }

    /** Post-order pass guaranteeing every parent is strictly older than its children. */
    private void enforceHeights(Node node) {
        if (node.isLeaf()) return;
        double maxChild = 0.0;
        for (Node child : node.getChildren()) {
            enforceHeights(child);
            maxChild = Math.max(maxChild, child.getHeight());
        }
        if (node.getHeight() <= maxChild) {
            node.setHeight(maxChild + MIN_GAP);
        }
    }

    // ------------------------------------------------------------------
    // Helpers

    private Node leaf(String name) {
        Node leaf = new Node();
        leaf.setID(name);
        leaf.setHeight(0.0);
        return leaf;
    }

    private Set<String> taxaOf(CalibrationNode node) {
        return taxaIds(node.taxa);
    }

    private static Set<String> taxaIds(TaxonSet ts) {
        Set<String> ids = new HashSet<>();
        for (Taxon t : ts.getTaxonSet()) ids.add(t.getID());
        return ids;
    }

    /** Number leaves from the taxon map and internal nodes sequentially, as {@code RandomTree} does. */
    private int numberNodes(Node node, int[] nextInternal, Map<String, Integer> taxonToNr) {
        if (node.isLeaf()) {
            Integer nr = taxonToNr.get(node.getID());
            if (nr == null) {
                throw new IllegalStateException("Unknown taxon in simulated tree: " + node.getID());
            }
            node.setNr(nr);
        } else {
            for (Node child : node.getChildren()) {
                numberNodes(child, nextInternal, taxonToNr);
            }
            node.setNr(taxaNames.size() + nextInternal[0]);
            nextInternal[0]++;
        }
        return nextInternal[0];
    }
}
