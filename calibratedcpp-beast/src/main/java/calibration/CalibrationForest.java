package calibration;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import calibrationprior.CalibrationCladePrior;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A forest of {@link CalibrationNode}s over a set of clades. The nesting structure is
 * derived purely from the clades' taxon sets; age bounds (carried by an optional
 * {@link CalibrationCladePrior} per node) are decoration that only some consumers read.
 *
 * <p>A forest may be built declaratively from XML via the {@code calibration} and/or
 * {@code taxonSet} inputs, programmatically via the constructors below, or parsed from a
 * Newick string by the {@link CalibrationForestParser} subclass (the analogue of
 * {@code TreeParser} to {@code Tree}).</p>
 *
 * @author Marcus Overwater
 */

@Description("A forest of calibration nodes.")
public class CalibrationForest extends BEASTObject {

    public Input<List<CalibrationCladePrior>> calibrationInput =
            new Input<>("calibration",
                    "Calibrated clades, each carrying a taxon set and soft lower/upper age bounds.",
                    new ArrayList<>());

    public Input<List<TaxonSet>> taxonSetInput =
            new Input<>("taxonSet",
                    "Taxon sets defining clades without age bounds (structure only). Clades whose "
                    + "taxa already appear in a 'calibration' entry are ignored here.",
                    new ArrayList<>());

    protected final List<CalibrationNode> allNodes = new ArrayList<>();
    protected final List<CalibrationNode> roots = new ArrayList<>();

    /** No-arg constructor for XML instantiation and subclassing (e.g. the Newick parser). */
    public CalibrationForest() {
    }

    /**
     * Builds a forest from a flat list of taxon sets. Nodes carry no calibration prior.
     */
    public CalibrationForest(List<TaxonSet> taxonSets) {
        for (TaxonSet ts : taxonSets) {
            allNodes.add(new CalibrationNode(ts));
        }
        buildInclusionForest();
    }

    /**
     * Builds a forest from a list of {@link CalibrationCladePrior} objects, wiring each
     * node's {@code prior} field so that {@link CalibrationNode#getCalibrationCladePrior()}
     * is available throughout.
     */
    public static CalibrationForest fromPriors(List<CalibrationCladePrior> priors) {
        List<TaxonSet> taxonSets = priors.stream()
                .map(CalibrationCladePrior::getTaxa)
                .collect(Collectors.toList());
        CalibrationForest forest = new CalibrationForest(taxonSets);
        // Wire prior references by matching taxon sets
        for (CalibrationNode node : forest.allNodes) {
            for (CalibrationCladePrior prior : priors) {
                if (taxonSetsEqual(node.taxa, prior.getTaxa())) {
                    node.prior = prior;
                    break;
                }
            }
        }
        return forest;
    }

    // --- Core forest structure builder ---
    private void buildInclusionForest() {
        roots.clear();

        // Precompute taxa sets
        Map<CalibrationNode, Set<String>> nodeTaxaMap = new HashMap<>();
        for (CalibrationNode node : allNodes) {
            nodeTaxaMap.put(node, taxaIds(node.taxa));
        }

        // --- Validate nesting/disjointness ---
        for (int i = 0; i < allNodes.size(); i++) {
            CalibrationNode a = allNodes.get(i);
            Set<String> aTaxa = nodeTaxaMap.get(a);

            for (int j = i + 1; j < allNodes.size(); j++) {
                CalibrationNode b = allNodes.get(j);
                Set<String> bTaxa = nodeTaxaMap.get(b);

                boolean aContainsB = aTaxa.containsAll(bTaxa);
                boolean bContainsA = bTaxa.containsAll(aTaxa);
                boolean disjoint = Collections.disjoint(aTaxa, bTaxa);

                if (aTaxa.equals(bTaxa)) {
                    throw new IllegalArgumentException("Duplicate calibration clade: " + a.taxa.getID());
                }

                if (!aContainsB && !bContainsA && !disjoint) {
                    throw new IllegalArgumentException(
                            "Calibration clades must be nested or disjoint: " +
                                    a.taxa.getID() + " and " + b.taxa.getID());
                }
            }
        }

        // --- Build parent-child relations ---
        for (CalibrationNode child : allNodes) {
            CalibrationNode bestParent = null;
            Set<String> childTaxa = nodeTaxaMap.get(child);

            for (CalibrationNode candidate : allNodes) {
                if (candidate == child) continue;

                Set<String> candidateTaxa = nodeTaxaMap.get(candidate);
                if (candidateTaxa.containsAll(childTaxa)) {
                    if (bestParent == null ||
                            (nodeTaxaMap.get(bestParent).containsAll(candidateTaxa)
                                    && !nodeTaxaMap.get(candidate).equals(childTaxa))) {
                        bestParent = candidate;
                    }
                }
            }

            if (bestParent == child) {
                throw new IllegalStateException("Self-parent detected for " + child);
            }

            child.parent = bestParent;
            child.isRoot = (bestParent == null);

            if (bestParent != null) {
                bestParent.children.add(child);
            } else {
                roots.add(child);
            }
        }

        // Detect cycles
        for (CalibrationNode node : allNodes) {
            checkForCycle(node, new HashSet<>());
        }
    }

    private void checkForCycle(CalibrationNode node, Set<CalibrationNode> visited) {
        if (node == null) return;
        if (!visited.add(node)) {
            throw new IllegalStateException("Cycle detected involving: " + node);
        }
        if (node.parent != null) {
            checkForCycle(node.parent, visited);
        }
    }

    // --- Lookup methods ---
    public CalibrationNode getCalibrationNodeFromTaxonSet(TaxonSet taxonSet) {
        Set<String> targetIds = taxaIds(taxonSet);
        for (CalibrationNode node : allNodes) {
            if (taxaIds(node.taxa).equals(targetIds)) return node;
        }
        return null;
    }

    // --- Accessors ---
    public List<CalibrationNode> getAllNodes() {
        return allNodes;
    }

    public List<CalibrationNode> getRoots() {
        return roots;
    }

    /** Taxon sets for every clade in the forest (bounded and unbounded), in forest order. */
    public List<TaxonSet> getTaxonSets() {
        List<TaxonSet> out = new ArrayList<>(allNodes.size());
        for (CalibrationNode node : allNodes) out.add(node.taxa);
        return out;
    }

    /** Calibration priors for the clades that carry age bounds, in forest order. */
    public List<CalibrationCladePrior> getCalibrationCladePriors() {
        List<CalibrationCladePrior> out = new ArrayList<>();
        for (CalibrationNode node : allNodes) {
            if (node.prior != null) out.add(node.prior);
        }
        return out;
    }

    @Override
    public void initAndValidate() {
        allNodes.clear();
        roots.clear();
        populateNodes();
        buildInclusionForest();
    }

    /**
     * Populates {@link #allNodes} from the declarative inputs. Subclasses (e.g. a Newick
     * parser) override this to source nodes from elsewhere; it is called before the inclusion
     * hierarchy is built. A clade supplied via {@code taxonSet} is skipped when a
     * {@code calibration} entry already covers the same taxa.
     */
    protected void populateNodes() {
        Set<Set<String>> seen = new HashSet<>();
        for (CalibrationCladePrior prior : calibrationInput.get()) {
            CalibrationNode node = new CalibrationNode(prior.getTaxa());
            node.prior = prior;
            allNodes.add(node);
            seen.add(taxaIds(prior.getTaxa()));
        }
        for (TaxonSet ts : taxonSetInput.get()) {
            if (seen.add(taxaIds(ts))) {
                allNodes.add(new CalibrationNode(ts));
            }
        }
    }

    // --- Helpers ---
    private static Set<String> taxaIds(TaxonSet ts) {
        return ts.getTaxonSet().stream().map(Taxon::getID).collect(Collectors.toSet());
    }

    private static boolean taxonSetsEqual(TaxonSet a, TaxonSet b) {
        return taxaIds(a).equals(taxaIds(b));
    }
}
