package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.model.BirthDeathModel;
import calibratedcpp.model.CoalescentPointProcessModel;
import calibration.CalibrationClade;
import calibration.CalibrationForest;
import calibration.CalibrationNode;

import java.util.*;
import java.util.function.IntUnaryOperator;

/**
 * Species tree distribution for a coalescent point process conditioned on the ages of the mrca of clade calibrations.
 *
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

    public Input<List<CalibrationClade>> calibrationsInput =
            new Input<>("calibrations", "Clade calibrations", new ArrayList<>());

    public Input<Boolean> conditionOnCalibrationsInput =
            new Input<>("conditionOnCalibrations", "Boolean if the likelihood is conditioned on the clade calibrations (Default: true). " +
                    "For large trees with many calibrations it is recommended to set this to false and use the exchange operator.", true);

    protected List<CalibrationClade> calibrations;
    protected CalibrationForest calibrationForest;
    protected List<CalibrationNode> calibrationNodes;
    protected boolean conditionOnCalibrations;
    protected boolean conditionOnRoot;

    protected TreeInterface tree;
    protected CoalescentPointProcessModel model;
    protected Double origin;
    protected double rootAge;
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
            calibrationForest = new CalibrationForest(calibrations);
            calibrationNodes = calibrationForest.getAllNodes();
        }

        rootAge = tree.getRoot().getHeight();

        if (!conditionOnRoot && origin < rootAge) {
            throw new RuntimeException("rootAge (" + rootAge + ") > origin (" + origin + "): Origin must be greater than root age.");
        }

        maxTime = conditionOnRoot ? rootAge : origin;
    }

    public double calculateUnConditionedTreeLogLikelihood(TreeInterface tree) {

        double logP = Math.log1p(-Math.exp(model.calculateLogCDF(maxTime)));
        int numTaxa = tree.getLeafNodeCount();

        for (Node node : tree.getInternalNodes()) {
            double age = node.getHeight();
            logP += model.calculateLogDensity(age);
        }
        logP += (numTaxa - 1) * Math.log(2.0) - logFactorial(numTaxa); // Ignore orientation and include labelling with factor of 2^{n-1}/n!
        return logP;
    }

    protected double computeCalibrationDensity(TreeInterface tree, CalibrationNode calibration) {
        // Sort children descending by height
        List<CalibrationNode> children = calibration.children;
        children.sort(Comparator.comparingDouble(c -> -c.getCommonAncestor(tree).getHeight())); // sorts children from oldest to youngest

        Node mrca = calibration.getCommonAncestor(tree);
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

        double[] logDiff = new double[numChildren];
        double[] weights = new double[numChildren];

        for (int i = 0; i < numChildren; i++) {
            CalibrationNode child = children.get(i);
            logChildDensities[i] = computeCalibrationDensity(tree, child);
            childCladeSizes[i] = child.taxa.getTaxonSet().size();
            childCDFs[i] = model.calculateLogCDF(child.getCommonAncestor(tree).getHeight());

            logDiff[i] = logDiffExp(logQ_t, childCDFs[i]); // log(Q(t)-Q(x_i))
            weights[i] = logDiff[i] - logQ_t;
        }

        int sumChildSizes = Arrays.stream(childCladeSizes).sum();

        // Compute permutation sum over root locations
        logDensity += computeExtendedRootSum(weights, cladeSize - sumChildSizes) - logFactorial(cladeSize);

        for (int i = 0; i < numChildren; i++) {
            logDensity += logFactorial(childCladeSizes[i]) + 2 * logDiff[i] + logChildDensities[i];
        }

        logDensity += logFactorial(cladeSize - sumChildSizes) + (cladeSize - sumChildSizes - 2 - numChildren) * logQ_t;

        return logDensity;
    }

    public double calculateMarginalLogDensityOfCalibrations(
            TreeInterface tree,
            CalibrationForest calibrationForest) {

        List<CalibrationNode> rootNodes = calibrationForest.getRoots();

        double logQ_t = model.calculateLogCDF(maxTime);

        int numRoots = rootNodes.size();
        double[] logRootDensities = new double[numRoots];
        int[] rootCladeSizes = new int[numRoots];
        double[] rootCDFs = new double[numRoots];

        double[] logDiff = new double[numRoots];
        double[] weights = new double[numRoots];

        for (int i = 0; i < numRoots; i++) {
            logRootDensities[i] = computeCalibrationDensity(tree, rootNodes.get(i));
            rootCladeSizes[i] = rootNodes.get(i).taxa.getTaxonSet().size();
            rootCDFs[i] = model.calculateLogCDF(rootNodes.get(i).getCommonAncestor(tree).getHeight());

            logDiff[i] = logDiffExp(logQ_t, rootCDFs[i]);
            weights[i] = logDiff[i] - logQ_t;
        }

        int numFreeLineages = tree.getLeafNodeCount() - Arrays.stream(rootCladeSizes).sum(); // number of free lineages

        double interactionSum;
        if (conditionOnRoot || numFreeLineages==0) {
            interactionSum = computeExtendedRootSum(weights, numFreeLineages) + model.calculateLogDensity(maxTime);
        } else {
            interactionSum = computeBellmanHeldKarpWithTruncatedESP(weights, numFreeLineages);
        }

        double density = interactionSum + (numFreeLineages - numRoots - 2) * logQ_t + logFactorial(numFreeLineages) - logFactorial(tree.getLeafNodeCount());

        for (int i = 0; i< numRoots; i++) {
            density += logFactorial(rootCladeSizes[i]) + logRootDensities[i] + 2 * logDiff[i];
        }

        return density;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        if (!conditionOnRoot && origin < rootAge) {
            return Double.NEGATIVE_INFINITY;
        }

        logP = 0.0;

        logP += calculateUnConditionedTreeLogLikelihood(tree);

        if (conditionOnCalibrations) {
            for (CalibrationNode c : calibrationNodes) {
                Set<String> leafIDs = new HashSet<>();
                Node mrca = c.getCommonAncestor(tree);

                collectLeafTaxa(mrca, leafIDs);
                if (!leafIDs.equals(c.taxa.getTaxaNames())) {
                    return Double.NEGATIVE_INFINITY; // clade is not monophyletic!
                }
                if (c.getCalibrationClade().providedAge) {
                    if (Math.abs(mrca.getHeight() - c.getCalibrationClade().getAge().getValue()) > 1e-4) {
                        return Double.NEGATIVE_INFINITY; // tmrca of clade is not the calibration age!
                    }
                }
            }
            logP -= calculateMarginalLogDensityOfCalibrations(tree, calibrationForest) + Math.log1p(-Math.exp(model.calculateLogCDF(maxTime)));
        }

        return logP;
    }

    private double logDiffExp(double a, double b) {
//        if (b > a) throw new IllegalArgumentException("logDiffExp: b must be <= a");
        if (b > a) return Double.NEGATIVE_INFINITY;
        if (a == b) return Double.NEGATIVE_INFINITY;
        return a + Math.log1p(-Math.exp(b - a));
    }

    private double logFactorial(double n) {
        double result = 0.0;
        if (n < 0) return Double.NEGATIVE_INFINITY;
        for (int i = 1; i <= n; i++) {
            result += Math.log(i);
        }
        return result;
    }

    private double logAdd(double logA, double logB) {
        if (Double.isInfinite(logA)) return logB;
        if (Double.isInfinite(logB)) return logA;
        double m = Math.max(logA, logB);
        return m + Math.log(Math.exp(logA - m) + Math.exp(logB - m));
    }

    public double computeBellmanHeldKarpWithTruncatedESP(double[] lx, int M) {
        int k = lx.length;
        if (k == 0) return 0.0;
        if (k >= 20) throw new IllegalArgumentException("k must be < 20 for this array-based implementation.");

        int W = M + 1;
        int fullMask = (1 << k) - 1;

        // compute total array size and guard integer overflow
        long totalStates = (long) (1 << k) * k * W;
        if (totalStates > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("DP array size too large: " + totalStates);
        }
        int totalSize = (int) totalStates;

        // flat dp array: index = ((mask * k) + j) * W + w
        double[] dp = new double[totalSize];
        Arrays.fill(dp, Double.NEGATIVE_INFINITY);

        // helper for indexing
        int strideMask = k * W; // number of entries per mask
        java.util.function.IntUnaryOperator idx = (mi) -> mi * strideMask; // base for mask

        // Precompute logs:
        // la[j] = log(a_j) = log(1/x_j) = -lx[j]
        double[] la = new double[k];
        for (int j = 0; j < k; j++) la[j] = -lx[j];

        // laij[i][j] = log(a_{i,j}) = log(1/max(x_i,x_j)) = -max(lx[i], lx[j])
        double[][] laij = new double[k][k];
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                laij[i][j] = -Math.max(lx[i], lx[j]);
            }
        }

        // Precompute log-factorials up to (M-1)
        int nForBinom = Math.max(0, M - 1);
        double[] logFact = new double[nForBinom + 1];
        logFact[0] = 0.0;
        for (int i = 1; i <= nForBinom; i++) logFact[i] = logFact[i - 1] + Math.log(i);

        // Precompute log binomials logBin[w] = log C(M-1, w-1) for w=0..M
        double[] logBin = new double[W];
        for (int w = 0; w <= M; w++) {
            if (M == 0 && w == 0) {
                // Special case: binom(-1, -1) = 1
                logBin[w] = 0.0; // log(1)
                continue;
            }
            else if (w == 0) {
                logBin[w] = Double.NEGATIVE_INFINITY; // normally C(M-1, -1)=0
                continue;
            }
            int r = w - 1;
            if (r < 0 || r > M - 1)
                logBin[w] = Double.NEGATIVE_INFINITY;
            else
                logBin[w] = logFact[M - 1] - logFact[r] - logFact[(M - 1) - r];
        }

        // Base cases: singletons
        for (int j = 0; j < k; j++) {
            int mask = 1 << j;
            int base = idx.applyAsInt(mask);
            // dp[((mask * k) + j) * W + w]
            int baseIdx = base + j * W;
            for (int w = 0; w < W; w++) dp[baseIdx + w] = Double.NEGATIVE_INFINITY;
            // P_{ {j}, j }(t) = a_j + t  => coeff[0] = a_j, coeff[1] = 1
            dp[baseIdx] = la[j];         // log(a_j)
            if (W > 1) dp[baseIdx + 1] = 0.0; // log(1)
        }

        // Iterate masks in increasing order
        int maxMask = 1 << k;
        for (int mask = 1; mask < maxMask; mask++) {
            int bits = Integer.bitCount(mask);
            if (bits == 1) continue; // already inited

            int baseMask = idx.applyAsInt(mask);

            // For each endpoint j in mask
            for (int j = 0; j < k; j++) {
                if ((mask & (1 << j)) == 0) continue; // j not in mask
                int baseIdx = baseMask + j * W;

                // initialize to Double.NEGATIVE_INFINITY
                for (int w = 0; w < W; w++) dp[baseIdx + w] = Double.NEGATIVE_INFINITY;

                int prevMask = mask ^ (1 << j);
                int basePrev = idx.applyAsInt(prevMask);

                // iterate over possible previous endpoint i in prevMask
                for (int i = 0; i < k; i++) {
                    if ((prevMask & (1 << i)) == 0) continue;
                    int prevIdx = basePrev + i * W;

                    // for each coefficient w
                    // dp[mask,j,w] += laij[i][j] + dp[prevMask,i,w]
                    // dp[mask,j,w] += dp[prevMask,i,w-1]   (for w>0)
                    for (int w = 0; w < W; w++) {
                        double cur = dp[baseIdx + w];

                        // term a_{i,j} * prev[w]  -> log: prev[w] + laij[i][j]
                        double prevW = dp[prevIdx + w];
                        if (!Double.isInfinite(prevW)) {
                            double v = prevW + laij[i][j];
                            cur = logAdd(cur, v);
                        }

                        // term prev[w-1] (if w>0)
                        if (w > 0) {
                            double prevWm1 = dp[prevIdx + (w - 1)];
                            if (!Double.isInfinite(prevWm1)) {
                                cur = logAdd(cur, prevWm1);
                            }
                        }

                        dp[baseIdx + w] = cur;
                    }
                }
            }
        }

        // Finalize: for each endpoint j, compute ~P_j[w] = a_j * P[full][j][w] + P[full][j][w-1]
        double total = 0.0;
        int baseFull = idx.applyAsInt(fullMask);
        for (int j = 0; j < k; j++) {
            int pIdx = baseFull + j * W;
            // build ~P_j in log-space, accumulate weighted sum
            for (int w = 0; w < W; w++) {
                double tLog = Double.NEGATIVE_INFINITY;
                double pw = dp[pIdx + w];
                if (!Double.isInfinite(pw)) tLog = logAdd(tLog, pw + la[j]); // a_j * P[w] => add logs

                if (w > 0) {
                    double pwm1 = dp[pIdx + (w - 1)];
                    if (!Double.isInfinite(pwm1)) tLog = logAdd(tLog, pwm1);
                }

                if (Double.isInfinite(tLog) || Double.isInfinite(logBin[w])) continue;
                // convert back and accumulate: exp(tLog + logBin[w])
                total = logAdd(total, tLog + logBin[w]);
            }
        }

        return total;
    }

    public double computeExtendedRootSum(double[] lx, int M) {
        int k = lx.length;
        if (k == 0) return Double.NEGATIVE_INFINITY; // log(0)
        if (k >= 20) throw new IllegalArgumentException("k must be < 20 for this array-based implementation.");

        int W = M + 1;
        int fullMask = (1 << k) - 1;

        // 1. Memory Management
        long totalStates = (long) (1 << k) * k * W;
        if (totalStates > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("DP array size too large.");
        }
        int totalSize = (int) totalStates;

        // dpP stores log(Sum(weights))
        double[] dpP = new double[totalSize];
        Arrays.fill(dpP, Double.NEGATIVE_INFINITY);

        // dpQ stores log(Sum(weights * cost))
        double[] dpQ = new double[totalSize];
        Arrays.fill(dpQ, Double.NEGATIVE_INFINITY);

        // Helper for flat indexing: index = ((mask * k) + j) * W + w
        int strideMask = k * W;
        IntUnaryOperator idx = (mi) -> mi * strideMask;

        // 2. Precomputations

        // lx contains log(x_j).
        // We need log(a_j) = log(1/x_j) = -log(x_j) = -lx[j].
        double[] la = new double[k];
        for (int j = 0; j < k; j++) la[j] = -lx[j];

        // laij[i][j] = log(1/min(x_i, x_j)) = -log(min(x_i, x_j)) = -min(lx[i], lx[j])
        double[][] laij = new double[k][k];

        // Cost 0: log(1/max(x_i, x_j)) = -log(max(x_i, x_j)) = -max(lx[i], lx[j])
        double[][] logCost0 = new double[k][k];

        // Cost 1: log(1/x_i + 1/x_j) = log( exp(-lx[i]) + exp(-lx[j]) )
        double[][] logCost1 = new double[k][k];

        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                // Weight term (min in log domain is just min)
                laij[i][j] = -Math.max(lx[i], lx[j]);

                // Cost terms
                logCost0[i][j] = -Math.min(lx[i], lx[j]);
                // Use logAdd for sum of exponentials
                logCost1[i][j] = logAdd(-lx[i], -lx[j]);
            }
        }

        // Log Factorials and Binomials
        int nForBinom = Math.max(0, M - 1);
        double[] logFact = new double[nForBinom + 1];
        logFact[0] = 0.0;
        for (int i = 1; i <= nForBinom; i++) logFact[i] = logFact[i - 1] + Math.log(i);

        double[] logBin = new double[W];
        for (int w = 0; w <= M; w++) {
            if (M == 0 && w == 0) { logBin[w] = 0.0; continue; }
            else if (w == 0) { logBin[w] = Double.NEGATIVE_INFINITY; continue; }
            int r = w - 1;
            if (r < 0 || r > M - 1) logBin[w] = Double.NEGATIVE_INFINITY;
            else logBin[w] = logFact[M - 1] - logFact[r] - logFact[(M - 1) - r];
        }

        // 3. Base Cases: Singletons {j}
        for (int j = 0; j < k; j++) {
            // If the weight is -Infinity (prob 0), skip or let it propagate as -Inf.

            int baseIdx = idx.applyAsInt(1 << j) + j * W;

            // Case z_0 = 0 (w=0): P = a_j, Q = P * (z_0 * ...) = 0 -> log(-inf)
            dpP[baseIdx] = la[j];
            dpQ[baseIdx] = Double.NEGATIVE_INFINITY;

            // Case z_0 = 1 (w=1): P = 1, Q = P * (a_j) = a_j
            if (W > 1) {
                dpP[baseIdx + 1] = 0.0; // log(1)
                dpQ[baseIdx + 1] = la[j];
            }
        }

        // 4. DP Iteration
        int maxMask = 1 << k;
        for (int mask = 1; mask < maxMask; mask++) {
            if (Integer.bitCount(mask) == 1) continue;

            int baseMask = idx.applyAsInt(mask);

            for (int j = 0; j < k; j++) {
                if ((mask & (1 << j)) == 0) continue;

                int baseIdx = baseMask + j * W;
                int prevMask = mask ^ (1 << j);
                int basePrev = idx.applyAsInt(prevMask);

                for (int i = 0; i < k; i++) {
                    if ((prevMask & (1 << i)) == 0) continue;
                    int prevIdx = basePrev + i * W;

                    // Compute transitions for all w
                    for (int w = 0; w < W; w++) {
                        double curP = dpP[baseIdx + w];
                        double curQ = dpQ[baseIdx + w];

                        // --- Case A: No Gap (z = 0), comes from same w ---
                        double prevP = dpP[prevIdx + w];
                        double prevQ = dpQ[prevIdx + w];

                        if (!Double.isInfinite(prevP)) {
                            // P_new = P_old * weight
                            double pTerm = prevP + laij[i][j];

                            // Q_new = Q_old * weight + P_old * weight * cost
                            double qTermPart1 = (!Double.isInfinite(prevQ)) ? prevQ + laij[i][j] : Double.NEGATIVE_INFINITY;
                            double qTermPart2 = pTerm + logCost0[i][j];

                            double qTerm = logAdd(qTermPart1, qTermPart2);

                            curP = logAdd(curP, pTerm);
                            curQ = logAdd(curQ, qTerm);
                        }

                        // --- Case B: Gap (z = 1), comes from w-1 ---
                        if (w > 0) {
                            double prevPm1 = dpP[prevIdx + (w - 1)];
                            double prevQm1 = dpQ[prevIdx + (w - 1)];

                            if (!Double.isInfinite(prevPm1)) {
                                // Q_new = Q_old + P_old * cost
                                // * 1
                                double qTermPart2 = prevPm1 + logCost1[i][j];
                                double qTerm = logAdd(prevQm1, qTermPart2);

                                curP = logAdd(curP, prevPm1);
                                curQ = logAdd(curQ, qTerm);
                            }
                        }

                        dpP[baseIdx + w] = curP;
                        dpQ[baseIdx + w] = curQ;
                    }
                }
            }
        }

        // 5. Final Summation
        double totalLogSum = Double.NEGATIVE_INFINITY;
        int baseFull = idx.applyAsInt(fullMask);

        for (int j = 0; j < k; j++) {
            int pIdx = baseFull + j * W;

            for (int w = 0; w < W; w++) {
                // Calculate ~P_j[w] and ~Q_j[w] incorporating the final edge z_k

                double pFinal = Double.NEGATIVE_INFINITY;
                double qFinal = Double.NEGATIVE_INFINITY;

                // Term z_k = 0: P accumulates a_j, Cost adds 0
                double pRaw = dpP[pIdx + w];
                double qRaw = dpQ[pIdx + w];
                if (!Double.isInfinite(pRaw)) {
                    pFinal = logAdd(pFinal, pRaw + la[j]);
                    // Q is scaled by a_j, no added cost
                    if (!Double.isInfinite(qRaw)) {
                        qFinal = logAdd(qFinal, qRaw + la[j]);
                    }
                }

                // Term z_k = 1: P accumulates 1, Cost adds a_j (requires w > 0)
                if (w > 0) {
                    double pRawM1 = dpP[pIdx + (w - 1)];
                    double qRawM1 = dpQ[pIdx + (w - 1)];
                    if (!Double.isInfinite(pRawM1)) {
                        pFinal = logAdd(pFinal, pRawM1);
                        // Q accumulates 1 + P * cost(a_j)
                        double termQ = (!Double.isInfinite(qRawM1)) ? qRawM1 : Double.NEGATIVE_INFINITY;
                        double termCost = pRawM1 + la[j];
                        qFinal = logAdd(qFinal, logAdd(termQ, termCost));
                    }
                }

                if (Double.isInfinite(pFinal) || Double.isInfinite(logBin[w])) continue;

                // Combine: Binom * ( Q_final + P_final * Constant )
                // Constant = M + k - w - 1
                int constantVal = M - w;

                double weightedPart; // log(P * C)
                boolean isNegative = false;

                if (constantVal == 0) {
                    weightedPart = Double.NEGATIVE_INFINITY;
                } else {
                    weightedPart = pFinal + Math.log(Math.abs(constantVal));
                    if (constantVal < 0) isNegative = true;
                }

                double termInside;
                if (!isNegative) {
                    termInside = logAdd(qFinal, weightedPart);
                } else {
                    // If P*C is negative in linear domain, we compute Q - |P*C|
                    termInside = logDiffExp(qFinal, weightedPart);
                }

                totalLogSum = logAdd(totalLogSum, logBin[w] + termInside);
            }
        }

        return totalLogSum;
    }

    // tree methods
    private void collectLeafTaxa(Node node, Set<String> leafIDs) {
        if (node.isLeaf()) {
            leafIDs.add(node.getID());
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                collectLeafTaxa(node.getChild(i), leafIDs);
            }
        }
    }

    public void updateModel() {
        model = cppModelInput.get();
        origin = (originInput.get() != null) ? originInput.get().getValue() : null;
        rootAge = tree.getRoot().getHeight();
        maxTime = (conditionOnRoot) ? rootAge : origin;
        calibrations = calibrationsInput.get();
    }

    @Override
    public boolean requiresRecalculation() {
        updateModel();
        return true;
    }

    @Override
    public void restore() {
        updateModel();
        super.restore();
    }

    public static void main(String[] args) {
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:2,B:2):1,C:3):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        BirthDeathModel birthDeath = new BirthDeathModel();

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
