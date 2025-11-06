package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.model.BirthDeathModel;
import calibratedcpp.model.CoalescentPointProcessModel;
import calibrationprior.CalibrationNode;

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

    public Input<List<TaxonSet>> calibrationsInput =
            new Input<>("calibrations", "Clade calibrations", new ArrayList<>());

    public Input<Boolean> conditionOnCalibrationsInput =
            new Input<>("conditionOnCalibrations", "Boolean if the likelihood is conditioned on the clade calibrations (Default: true). " +
                    "For large trees with many calibrations it is recommended to set this to false and use the exchange operator.", true);

    protected TreeInterface tree;

    protected CoalescentPointProcessModel model;
    protected List<TaxonSet> calibrations;
    protected List<CalibrationNode> calibrationForest;
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
            calibrationForest = buildCalibrationForest(tree, calibrations);
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

    protected double computeCalibrationDensity(TreeInterface tree, CalibrationNode calibration) {
        updateModel(tree);

        // Sort children descending by height
        List<CalibrationNode> children = calibration.children;
        children.sort(Comparator.comparingDouble(c -> -c.mrca.getHeight()));

        Node mrca = calibration.mrca;
        double cladeHeight = mrca.getHeight();
        int cladeSize = calibration.taxa.getTaxonSet().size();

        double logQ_t = model.calculateLogCDF(cladeHeight);
        double logDensity = model.calculateLogDensity(cladeHeight);

        // Base case: leaf calibration
        if (children.isEmpty()) {
            logDensity += (cladeSize - 2) * logQ_t + Math.log(cladeSize - 1);
            return logDensity;
        }

        // Recursive case: internal calibration
        int numChildren = children.size();
        double[] logChildDensities = new double[numChildren];
        int[] childCladeSizes = new int[numChildren];
        double[] childCDFs = new double[numChildren];

        for (int i = 0; i < numChildren; i++) {
            CalibrationNode child = children.get(i);
            logChildDensities[i] = computeCalibrationDensity(tree, child);
            childCladeSizes[i] = child.taxa.getTaxonSet().size();
            childCDFs[i] = model.calculateLogCDF(child.mrca.getHeight());
        }

        double[] logDiff = new double[numChildren];
        for (int i = 0; i < numChildren; i++) {
            logDiff[i] = logDiffExp(logQ_t, childCDFs[i]);
        }

        int sumChildSizes = Arrays.stream(childCladeSizes).sum();
        int numRootLocations = cladeSize - sumChildSizes + numChildren - 1;

        // Compute permutation sum over root locations
        List<Double> perRootLogs = new ArrayList<>();
        for (int rootLoc = 1; rootLoc <= numRootLocations; rootLoc++) {
            perRootLogs.add(calculateLogSumOfPermutationsWithRoot(
                    numChildren,
                    sumChildSizes,
                    cladeSize,
                    logQ_t,
                    logDiff,
                    rootLoc
            ));
        }

        double logPermutationSum = logSumExp(perRootLogs);
        double childrenSum = Arrays.stream(logChildDensities).sum();

        return logDensity + childrenSum + logPermutationSum;
    }

    public double calculateLogMarginalDensityOfCalibrations(
            TreeInterface tree,
            List<CalibrationNode> calibrations) {

        // Step 1: Update model and get base quantities
        updateModel(tree);
        double logQt = model.calculateLogCDF(maxTime);
        double marginalDensity = Math.log1p(-Math.exp(logQt)); // terminating node age > maxTime

        int numTaxa = tree.getLeafNodeCount();

        // Step 2: Collect and sort root calibrations by clade height (descending)
        List<CalibrationNode> roots = calibrations.stream()
                .filter(c -> c.isRoot)
                .sorted(Comparator.comparingDouble(c -> -c.mrca.getHeight()))
                .toList();

        // Step 3: Compute density for each root recursively
        int numRoots = roots.size();
        double[] logQti = new double[numRoots];
        double[] logDiff = new double[numRoots];
        int[] cladeSizes = new int[numRoots];
        int sumCladeSizes = 0;

        for (int i = 0; i < numRoots; i++) {
            CalibrationNode root = roots.get(i);
            double logDensity = computeCalibrationDensity(tree, root);
            marginalDensity += logDensity;

            double height = root.mrca.getHeight();
            logQti[i] = model.calculateLogCDF(height);
            logDiff[i] = logDiffExp(logQt, logQti[i]);

            cladeSizes[i] = root.taxa.getTaxonCount();
            sumCladeSizes += cladeSizes[i];
        }

        // Step 4: Add permutation terms depending on conditioning
        if (!conditionOnRoot) {
            marginalDensity += calculateLogSumOfPermutations(
                    numRoots, sumCladeSizes, numTaxa, logQt, logDiff);
        } else {
            int numRootLocations = numTaxa - sumCladeSizes + numRoots - 1;
            marginalDensity += Math.log1p(-Math.exp(logQt));
            for (int rootLoc = 1; rootLoc <= numRootLocations; rootLoc++) {
                marginalDensity += calculateLogSumOfPermutationsWithRoot(
                        numRoots, sumCladeSizes, numTaxa, logQt, logDiff, rootLoc);
            }
        }

        return marginalDensity;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        updateModel(tree);
        if (!conditionOnRoot && origin < rootAge) {
            return Double.NEGATIVE_INFINITY;
        }

        logP = 0.0;
        for (CalibrationNode c : calibrationForest) {
            Set<String> leafIDs = new HashSet<>();
            Node beastNode = c.mrca;

            collectLeafTaxa(beastNode, leafIDs);
            if (!leafIDs.equals(c.taxa.getTaxaNames())){
                return Double.NEGATIVE_INFINITY; // clade is not monophyletic!
            }
        }
        logP += calculateUnConditionedTreeLogLikelihood(tree);

        if (conditionOnCalibrations) {
            logP -= calculateLogMarginalDensityOfCalibrations(tree, calibrationForest);
        }

        return logP;
    }

    private double calculateLogSumOfPermutationsCommon(
            int numCalibrations,
            int sumCladeSizes,
            int numTaxa,
            double logQ_t,
            double[] logDiff,
            boolean hasRoot,
            int locationOfRoot) {

        int m = numTaxa - sumCladeSizes;
        int totalElements = m + numCalibrations;

        List<List<Integer>> permutations = new ArrayList<>();
        generatePermutations(totalElements, numCalibrations, new ArrayList<>(),
                new boolean[totalElements], permutations);

        List<Double> logTerms = new ArrayList<>();

        for (List<Integer> perm : permutations) {
            double logTerm = 0.0;
            int sum_s = 0;

            for (int i = 0; i < numCalibrations; i++) {
                int ell_i = perm.get(i);

                // count adjacency to previously chosen elements
                int countAdj = 0;
                for (int j = 0; j < i; j++) {
                    int ell_j = perm.get(j);
                    if (ell_j == ell_i - 1 || ell_j == ell_i + 1) {
                        countAdj++;
                    }
                }

                int s_i = 0;
                if (ell_i == 0 || ell_i == totalElements - 1) s_i++;
                s_i += countAdj;

                // Root adjacency logic
                if (hasRoot && (ell_i == locationOfRoot || ell_i == locationOfRoot - 1)) {
                    s_i++;
                }

                if (hasRoot) {
                    s_i = Math.min(s_i, 2);
                }

                sum_s += s_i;
                logTerm += (2 - s_i) * logDiff[i];
            }

            // Base term
            int rootAdjustment = hasRoot ? -1 : 0;
            logTerm += (m - numCalibrations - 1 + sum_s + rootAdjustment) * logQ_t;

            logTerms.add(logTerm);
        }

        return logSumExp(logTerms);
    }

    private double calculateLogSumOfPermutations(
            int numCalibrations, int sumCladeSizes, int numTaxa,
            double logQ_t, double[] logDiff) {

        return calculateLogSumOfPermutationsCommon(
                numCalibrations, sumCladeSizes, numTaxa,
                logQ_t, logDiff,
                false, -1);
    }

    private double calculateLogSumOfPermutationsWithRoot(
            int numCalibrations, int sumCladeSizes, int numTaxa,
            double logQ_t, double[] logDiff, int locationOfRoot) {

        return calculateLogSumOfPermutationsCommon(
                numCalibrations, sumCladeSizes, numTaxa,
                logQ_t, logDiff,
                true, locationOfRoot);
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

    // Log tricks
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

    // tree methods
    private boolean isDescendant(Node node, Node ancestor) {
        while (node != null) {
            if (node == ancestor) return true;
            node = node.getParent();
        }
        return false;
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

    // build calibration forest
    protected List<CalibrationNode> buildCalibrationForest(TreeInterface tree, List<TaxonSet> calibrations){
        List<CalibrationNode> calibrationNodes = new ArrayList<>();
        for (TaxonSet calibration : calibrations) {
            CalibrationNode calibrationNode = new CalibrationNode(tree, calibration);
            calibrationNodes.add(calibrationNode);
        }

        for (CalibrationNode child : calibrationNodes) {
            CalibrationNode bestParent = null;
            for (CalibrationNode cand : calibrationNodes) {
                if (cand == child) continue;
                if (isDescendant(child.mrca, cand.mrca)) {
                    // choose the closest ancestor, not the highest
                    if (bestParent == null || isDescendant(cand.mrca, bestParent.mrca)) {
                        bestParent = cand;
                    }
                }
            }
            child.parent = bestParent;
            if (bestParent != null) {
                bestParent.children.add(child);
                child.isRoot = false;
            } else {
                child.isRoot = true;
            }
        }
        return calibrationNodes;
    }

    public void updateModel(TreeInterface tree) {
        model = cppModelInput.get();
        origin = (originInput.get() != null) ? originInput.get().getValue() : null;
        rootAge = tree.getRoot().getHeight();
        maxTime = (conditionOnRoot) ? rootAge : origin;

        calibrations = new ArrayList<>(calibrationsInput.get());
        conditionOnCalibrations = (!calibrations.isEmpty()) ? conditionOnCalibrationsInput.get() : false;

        if (conditionOnCalibrations) {
            calibrationForest = buildCalibrationForest(tree, calibrations);
        }
    }

    @Override
    public boolean requiresRecalculation() {
        buildCalibrationForest(treeInput.get(), calibrations);
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
