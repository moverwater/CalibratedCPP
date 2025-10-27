package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.model.BirthDeathModel;
import calibratedcpp.model.CoalescentPointProcessModel;

import java.util.*;

/**
 * @author Marcus Overwater
 */

@Description("A general class of birth-death processes with incomplete extant sampling and conditioning on clade calibrations")
public class CalibratedCoalescentPointProcess extends SpeciesTreeDistribution {
    public Input<RealParameter> originInput =
            new Input<>("origin", "Age of the origin (time of process start)", (RealParameter) null);

    public Input<Boolean> conditionOnRootInput =
            new Input<>("conditionOnRoot", "Whether the model is conditioned on the root age (default: false)", false);

    public Input<CoalescentPointProcessModel> cppModelInput =
            new Input<>("treeModel", "The tree model", (CoalescentPointProcessModel) null);

    public Input<List<CalibrationPoint>> calibrationsInput =
            new Input<>("calibrations", "Clade calibrations", new ArrayList<>());

    public Input<Boolean> conditionOnCalibrationsInput =
            new Input<>("conditionOnCalibrations", "Boolean if the likelihood is conditioned on the clade calibrations (Default: true). " +
                    "For large trees with many calibrations it is recommended to set this to false and use the exchange operator.", true);

    protected TreeInterface tree;

    protected CoalescentPointProcessModel model;
    protected List<CalibrationPoint> calibrations;
    protected Map<CalibrationPoint, List<CalibrationPoint>> calibrationGraph;
    protected boolean conditionOnCalibrations;

    protected Double origin;
    protected Double rootAge;

    protected boolean conditionOnRoot;
    protected double maxTime;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        RealParameter originParam = originInput.get();
        origin = (originParam != null) ? originParam.getValue() : null;
        conditionOnRoot = conditionOnRootInput.get();

        if (origin == null && !conditionOnRoot) {
            throw new IllegalArgumentException("You must either provide an origin age or set conditionOnRoot=true.");
        }

        tree = treeInput.get();
        model = cppModelInput.get();
        calibrations = new ArrayList<>(calibrationsInput.get());
        conditionOnCalibrations = (!calibrations.isEmpty()) ? conditionOnCalibrationsInput.get() : false;

        if (conditionOnCalibrations) {
            calibrations = postOrderTopologicalSort(tree, calibrations);
            calibrationGraph = buildNestingDAG(calibrations, tree); // <-- THIS LINE IS MISSING
        } else {
            calibrationGraph = new HashMap<>();
        }

        rootAge = tree.getRoot().getHeight();

        if (!conditionOnRoot && origin < rootAge) {
            throw new RuntimeException("rootAge (" + rootAge + ") > origin (" + origin + "): Origin must be greater than root age.");
        }

        maxTime = conditionOnRoot ? rootAge : origin;
    }

    public double calculateUnConditionedTreeLogLikelihood(TreeInterface tree) {
        updateModel(tree);

        double logP = Math.log1p(-Math.exp(model.calculateLogCDF(maxTime)));

        if (conditionOnRoot) {
            logP += logP - model.calculateLogDensity(rootAge);
        }

        for (Node node : tree.getInternalNodes()) {
            double age = node.getHeight();
            logP += model.calculateLogDensity(age);
        }
        return logP;
    }

    public double calculateLogMarginalDensityOfCalibrations(TreeInterface tree,
                                                            List<CalibrationPoint> calibrations,
                                                            Map<CalibrationPoint, List<CalibrationPoint>> calibrationGraph) {

        updateModel(tree);
        double logQt = model.calculateLogCDF(maxTime);
        double marginalDensity = logDiffExp(0.0, logQt); // terminating node age greater than maxTime

        int numTaxa = tree.getLeafNodeCount();

        // Step 1: Find all CalibrationPoints that appear as children
        Set<CalibrationPoint> children = new HashSet<>();
        for (List<CalibrationPoint> childList : calibrationGraph.values()) {
            children.addAll(childList);
        }

        // Step 2: Roots = all calibrations that are not children
        Set<CalibrationPoint> roots = new HashSet<>(calibrations);
        roots.removeAll(children); // Now only root calibrations remain

        // Step 2b: Sort roots by clade height (oldest first)
        List<CalibrationPoint> sortedRoots = new ArrayList<>(roots);
        sortedRoots.sort((c1, c2) -> {
            double h1 = getMRCA(tree, c1.taxa().asStringList()).getHeight();
            double h2 = getMRCA(tree, c2.taxa().asStringList()).getHeight();
            return Double.compare(h2, h1); // descending: oldest first
        });

        // Step 3: Compute log density for each root (recursively handles children)
        int numMaxCalibrations = roots.size();

        double[] calibrationAges = new double[numMaxCalibrations];
        int[] cladeSizes = new int[numMaxCalibrations];
        int sumCladeSizes = 0;

        double[] logQti = new double[numMaxCalibrations];
        double[] logDiff = new double[numMaxCalibrations];

        int index = 0;
        for (CalibrationPoint root : sortedRoots) {
            marginalDensity += calculateLogDensityOfSingleCalibration(tree, root, calibrationGraph);

            calibrationAges[index] = getMRCA(tree, root.taxa().asStringList()).getHeight();
            logQti[index] = model.calculateLogCDF(calibrationAges[index]);
            logDiff[index] = logDiffExp(logQt, logQti[index]);

            cladeSizes[index] = root.taxa().getTaxonCount();
            sumCladeSizes += cladeSizes[index];

            index++;
        }

        if (!conditionOnRoot) {
            marginalDensity += calculateLogSumOfPermutations(numMaxCalibrations, sumCladeSizes, numTaxa, logQt, logDiff);
        } else {
            int numRootLocations = numTaxa - sumCladeSizes + numMaxCalibrations - 1;
            marginalDensity += logDiffExp(0.0, logQt);
            for (int rootLocation = 1; rootLocation <= numRootLocations; rootLocation++) {
                marginalDensity += calculateLogSumOfPermutationsWithRoot(numMaxCalibrations, sumCladeSizes, numTaxa, logQt, logDiff, rootLocation);
            }
        }

        return marginalDensity;
    }

    public double calculateLogDensityOfSingleCalibration(TreeInterface tree,
                                                         CalibrationPoint calibration,
                                                         Map<CalibrationPoint, List<CalibrationPoint>> calibrationGraph) {
        updateModel(tree);
        List<CalibrationPoint> children = calibrationGraph.getOrDefault(calibration, new ArrayList<>());
        children.sort((c1, c2) -> {
            double h1 = getMRCA(tree, c1.taxa().asStringList()).getHeight();
            double h2 = getMRCA(tree, c2.taxa().asStringList()).getHeight();
            return Double.compare(h2, h1); // descending: oldest (largest height) first
        });

        Node mrca = getMRCA(tree, calibration.taxa().asStringList());
        double cladeHeight = mrca.getHeight();
        int cladeSize = calibration.taxa().getTaxonSet().size();

        double logQ_t = model.calculateLogCDF(cladeHeight);
        double logDensity = model.calculateLogDensity(cladeHeight);

        if (children.isEmpty()) {
            // Leaf calibration
            logDensity += (cladeSize - 2) * logQ_t + Math.log(cladeSize - 1);
            return logDensity;
        }

        // Internal calibration with nested calibrations
        int numChildren = children.size();
        double[] logChildDensities = new double[numChildren];
        int[] childCladeSizes = new int[numChildren];
        double[] childCDFs = new double[numChildren];

        for (int i = 0; i < numChildren; i++) {
            CalibrationPoint child = children.get(i);
            logChildDensities[i] = calculateLogDensityOfSingleCalibration(tree, child, calibrationGraph);
            childCladeSizes[i] = child.taxa().getTaxonSet().size();
            childCDFs[i] = model.calculateLogCDF(getMRCA(tree, child.taxa().asStringList()).getHeight());
        }

        double[] logDiff = new double[numChildren];
        for (int i = 0; i < numChildren; i++) {
            logDiff[i] = logDiffExp(logQ_t, childCDFs[i]);
        }

        double logPermutationSum;
        int numRootLocations = cladeSize - sum(childCladeSizes) + numChildren - 1;
        List<Double> perRootLogs = new ArrayList<>();
        for (int rootLocation = 1; rootLocation < numRootLocations + 1; rootLocation++) {
            double val = calculateLogSumOfPermutationsWithRoot(
                    numChildren,
                    sum(childCladeSizes),
                    cladeSize,
                    logQ_t,
                    logDiff,
                    rootLocation
            );
            perRootLogs.add(val);
        }

        logPermutationSum = logSumExp(perRootLogs);

        double childrenTerm = 0.0;
        for (double logD : logChildDensities) {
            childrenTerm += logD;
        }

        logDensity += childrenTerm + logPermutationSum;
        return logDensity;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        updateModel(tree);
        if (!conditionOnRoot && origin < rootAge) {
            return Double.NEGATIVE_INFINITY;
        }

        logP = 0.0;
        for (CalibrationPoint c : calibrations) {
            Node mrca = getMRCA(tree, c.taxa().asStringList());
            Set<String> mrcaLeafIDs = new HashSet<>();
            collectLeafTaxa(mrca, mrcaLeafIDs);
            Set<String> calibrationTaxa = new HashSet<>(c.taxa().asStringList());

            if (!mrcaLeafIDs.equals(calibrationTaxa)) {
                return Double.NEGATIVE_INFINITY; // Calibration clade is not monophyletic
            }
        }
        logP += calculateUnConditionedTreeLogLikelihood(tree);

        if (conditionOnCalibrations) {
            logP -= calculateLogMarginalDensityOfCalibrations(tree, calibrations, calibrationGraph);
        }

        return logP;
    }

    private void collectLeafTaxa(Node node, Set<String> leafIDs) {
        if (node.isLeaf()) {
            leafIDs.add(node.getID());
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                collectLeafTaxa(node.getChild(i), leafIDs);
            }
        }
    }

    private double calculateLogSumOfPermutations(int numCalibrations, int sumCladeSizes, int numTaxa, double logQ_t, double[] logDiff) {

        int m = numTaxa - sumCladeSizes;

        int totalElements = m + numCalibrations;

        List<List<Integer>> permutations = new ArrayList<>();
        generatePermutations(totalElements, numCalibrations, new ArrayList<>(), new boolean[totalElements], permutations);

        List<Double> logTerms = new ArrayList<>();

        for (List<Integer> perm : permutations) {
            double logTerm = 0;
            int sum_s = 0;
            int[] s = new int[numCalibrations];

            for (int i = 0; i < numCalibrations; i++) {
                int ell_i = perm.get(i);
                int countAdj = 0;
                for (int j = 0; j < i; j++) {
                    int ell_j = perm.get(j);
                    if (ell_j == ell_i - 1 || ell_j == ell_i + 1) {
                        countAdj++;
                    }
                }
                s[i] = (ell_i == 0 ? 1 : 0) + (ell_i == totalElements - 1 ? 1 : 0) + countAdj;
                sum_s += s[i];

                logTerm += (2 - s[i]) * logDiff[i]; // Compute log term for this permutation
            }

            logTerm += (m - numCalibrations - 1 + sum_s) * logQ_t;

            logTerms.add(logTerm);
        }
        // Use logSumExp to sum all terms in log space
        return logSumExp(logTerms);
    }

    private double calculateLogSumOfPermutationsWithRoot(int numCalibrations, int sumCladeSizes, int numTaxa, double logQ_t, double[] logDiff, int locationOfRoot) {

        int m = numTaxa - sumCladeSizes;

        int totalElements = m + numCalibrations;

        List<List<Integer>> permutations = new ArrayList<>();
        generatePermutations(totalElements, numCalibrations, new ArrayList<>(), new boolean[totalElements], permutations);

        List<Double> logTerms = new ArrayList<>();

        for (List<Integer> perm : permutations) {
            double logTerm = 0;
            int sum_s = 0;
            int[] s = new int[numCalibrations];

            for (int i = 0; i < numCalibrations; i++) {
                int ell_i = perm.get(i);
                int countAdj = 0;
                for (int j = 0; j < i; j++) {
                    int ell_j = perm.get(j);
                    if (ell_j == ell_i - 1 || ell_j == ell_i + 1) {
                        countAdj++;
                    }
                }
                s[i] = (ell_i == 0 ? 1 : 0) + (ell_i == totalElements - 1 ? 1 : 0) + countAdj;
                // Extra adjacency if next to root location (but not at the ends)
                if (ell_i == locationOfRoot || ell_i == locationOfRoot - 1) {
                    s[i]++;
                }

                s[i] = Math.min(s[i], 2);
                sum_s += s[i];

                logTerm += (2 - s[i]) * logDiff[i]; // Compute log term for this permutation
            }

            logTerm += (m - numCalibrations - 1 + sum_s - 1) * logQ_t;

            logTerms.add(logTerm);
        }
        // Use logSumExp to sum all terms in log space
        return logSumExp(logTerms);
    }

    private Node getMRCA(TreeInterface tree, List<String> taxonIDs) {
        // Get leaf nodes for the given taxa
        List<Node> nodes = new ArrayList<>();
        for (String id : taxonIDs) {
            int index = tree.getTaxonset().getTaxonIndex(id);
            nodes.add(tree.getNode(index));
        }

        return findMRCA(nodes);
    }

    private Node findMRCA(List<Node> nodes) {
        if (nodes == null || nodes.isEmpty()) return null;
        if (nodes.size() == 1) return nodes.get(0);

        Node mrca = nodes.get(0);

        for (int i = 1; i < nodes.size(); i++) {
            mrca = findMRCA(mrca, nodes.get(i));
        }

        return mrca;
    }

    private Node findMRCA(Node a, Node b) {
        // Collect ancestors of a
        List<Node> ancestorsA = new ArrayList<>();
        while (a != null) {
            ancestorsA.add(a);
            a = a.getParent();
        }

        // Walk b's ancestry and return first shared node
        while (b != null) {
            for (Node ancestor : ancestorsA) {
                if (b == ancestor) { // use reference equality!
                    return b;
                }
            }
            b = b.getParent();
        }

        return null; // Should never happen if both are in same tree
    }

    private void generatePermutations(int n, int k, List<Integer> current, boolean[] used, List<List<Integer>> result) {
        if (current.size() == k) {
            result.add(new ArrayList<>(current));
            return;
        }

        for (int i = 0; i < n; i++) {
            if (!used[i]) {
                used[i] = true;
                current.add(i);
                generatePermutations(n, k, current, used, result);
                current.remove(current.size() - 1);
                used[i] = false;
            }
        }
    }

    private double logSumExp(List<Double> logValues) {
        if (logValues.isEmpty()) return Double.NEGATIVE_INFINITY;
        double max = Collections.max(logValues);
        double sum = 0.0;
        for (double val : logValues) sum += Math.exp(val - max);
        return max + Math.log(sum);
    }

    private double logDiffExp(double a, double b) {
//        if (b > a) throw new IllegalArgumentException("logDiffExp: b must be <= a");
        if (b > a) return Double.NEGATIVE_INFINITY;
        if (a == b) return Double.NEGATIVE_INFINITY;
        return a + Math.log1p(-Math.exp(b - a));
    }

    protected boolean isNested(CalibrationPoint inner, CalibrationPoint outer, TreeInterface tree) {
        Node innerMRCA = getMRCA(tree, inner.taxa().asStringList());
        Node outerMRCA = getMRCA(tree, outer.taxa().asStringList());
        return isDescendant(innerMRCA, outerMRCA);
    }

    private boolean isDescendant(Node node, Node ancestor) {
        while (node != null) {
            if (node == ancestor) return true;
            node = node.getParent();
        }
        return false;
    }

    protected Map<CalibrationPoint, List<CalibrationPoint>> buildNestingDAG(List<CalibrationPoint> calibrations, TreeInterface tree) {
        Map<CalibrationPoint, List<CalibrationPoint>> graph = new HashMap<>();

        // Initialize adjacency list
        for (CalibrationPoint c : calibrations) {
            graph.put(c, new ArrayList<>());
        }

        for (CalibrationPoint parent : calibrations) {
            for (CalibrationPoint child : calibrations) {
                if (parent == child) continue;
                if (isNested(child, parent, tree)) {
                    // Check if there is no intermediate calibration between parent and child
                    boolean hasIntermediate = false;
                    for (CalibrationPoint other : calibrations) {
                        if (other != parent && other != child &&
                                isNested(child, other, tree) && isNested(other, parent, tree)) {
                            hasIntermediate = true;
                            break;
                        }
                    }
                    if (!hasIntermediate) {
                        graph.get(parent).add(child);
                    }
                }
            }
        }

        // Sort children lists by MRCA height: oldest (largest height) first
        Map<CalibrationPoint, Double> cachedHeights = new HashMap<>();
        for (CalibrationPoint c : calibrations) {
            Node mrca = getMRCA(tree, c.taxa().asStringList());
            cachedHeights.put(c, (mrca == null ? Double.NEGATIVE_INFINITY : mrca.getHeight()));
        }

        for (List<CalibrationPoint> children : graph.values()) {
            children.sort((c1, c2) -> {
                double h1 = cachedHeights.getOrDefault(c1, Double.NEGATIVE_INFINITY);
                double h2 = cachedHeights.getOrDefault(c2, Double.NEGATIVE_INFINITY);
                return Double.compare(h2, h1); // descending
            });
        }

        return graph;
    }

    protected List<CalibrationPoint> postOrderTopologicalSort(TreeInterface tree, List<CalibrationPoint> calibrations) {
        Map<CalibrationPoint, List<CalibrationPoint>> graph = buildNestingDAG(calibrations, tree);
        List<CalibrationPoint> sorted = new ArrayList<>();
        Set<CalibrationPoint> visited = new HashSet<>();

        for (CalibrationPoint c : calibrations) {
            postOrderDFS(c, graph, visited, sorted, tree);
        }

        return sorted;
    }

    private void postOrderDFS(CalibrationPoint current, Map<CalibrationPoint, List<CalibrationPoint>> graph,
                              Set<CalibrationPoint> visited, List<CalibrationPoint> result, TreeInterface tree) {
        if (visited.contains(current)) return;
        visited.add(current);

        List<CalibrationPoint> children = graph.getOrDefault(current, new ArrayList<>());
        children.sort(Comparator.comparingDouble(
                (CalibrationPoint c) -> getMRCA(tree, c.taxa().asStringList()).getHeight()
        ).reversed()); // Oldest to youngest

        for (CalibrationPoint child : children) {
            postOrderDFS(child, graph, visited, result, tree);
        }

        result.add(current); // Post-order: visit after children
    }

    private int sum(int[] values) {
        int output = 0;
        for (int val : values) output += val;
        return output;
    }

    public void updateModel(TreeInterface tree) {
        model = cppModelInput.get();
        origin = (originInput.get() != null) ? originInput.get().getValue() : null;
        rootAge = tree.getRoot().getHeight();
        maxTime = (conditionOnRoot) ? rootAge : origin;

        calibrations = new ArrayList<>(calibrationsInput.get());
        conditionOnCalibrations = (!calibrations.isEmpty()) ? conditionOnCalibrationsInput.get() : false;

        if (conditionOnCalibrations) {
            calibrations = postOrderTopologicalSort(tree, calibrations);
            calibrationGraph = buildNestingDAG(calibrations, tree);
        } else {
            calibrationGraph = new HashMap<>();
        }
    }

    @Override
    public boolean requiresRecalculation() {
        super.requiresRecalculation();
        return true;
    }

    public static void main(String[] args) {
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        calibratedcpp.model.BirthDeathModel birthDeath = new BirthDeathModel();

        birthDeath.initByName("birthRate", new RealParameter("3.0"),
                "deathRate", new RealParameter("2.0"),
                "rho", new RealParameter("0.1")
        );

        boolean b;
        b = birthDeath.birthRateInput.get().getValue() - birthDeath.deathRateInput.get().getValue() < 1e-10;

        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "treeModel", birthDeath,
                "origin", new RealParameter("4.0")
        );
        System.out.println("tree = " + tree);
        System.out.println("isCritical = " + b);
        System.out.println("logP = " + cpp.calculateLogP());
    }
}
