package calibratedcpp.lphy.prior;

import calibratedcpp.lphy.util.TruncatedLogNormal;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;

import static calibratedcpp.lphy.prior.ConditionedPriorUtils.*;
import static calibratedcpp.lphy.tree.CPPUtils.*;

public class ConditionedMRCAPrior implements GenerativeDistribution<Calibration[]> {
    Value<String[][]> calibrationTaxa;
    Value<Number[]> upperBounds;
    Value<Number[]> lowerBounds;
    Value<Number> coverage;
    Value<Boolean> rootFlag;

    public final String calibrationTaxaName = "calibrationTaxa";
    public final String rootFlagName = "rootFlag";
    public final String upperBoundName = "upperBounds";
    public final String lowerBoundName = "lowerBounds";
    public final String coverageName = "p";

    public ConditionedMRCAPrior(
            @ParameterInfo(name = calibrationTaxaName, description = "the array of calibration names, the upper and lower bounds should in this order, root calibration should be put at the first") Value<String[][]> calibrationTaxa,
            @ParameterInfo(name = rootFlagName, description = "if root is also calibrated, default no", optional = true) Value<Boolean> rootFlag,
            @ParameterInfo(name = upperBoundName, description = "the array of the upper bounds of the corresponding calibration, if root is calibrated then put it at the first value") Value<Number[]> upperBounds,
            @ParameterInfo(name = lowerBoundName, description = "the array of the lower bounds of the corresponding calibration, if root is calibrated then put it at the first value") Value<Number[]> lowerBounds,
            @ParameterInfo(name = coverageName, description = "the confidential level that the amount of probability mass expected to be retained after truncation, default 0.9", optional = true) Value<Number> coverage
    ) {
        // check illegal arguments
        if (calibrationTaxa == null || calibrationTaxa.value() == null || calibrationTaxa.value().length < 1)
            throw new IllegalArgumentException("The calibrations must be equal or more than one");
        if (upperBounds == null || upperBounds.value() == null || upperBounds.value().length < 1)
            throw new IllegalArgumentException("The upper bounds must be equal or more than one");
        if (lowerBounds == null || lowerBounds.value() == null || lowerBounds.value().length < 1)
            throw new IllegalArgumentException("The lower bounds must be equal or more than one");

        if (calibrationTaxa.value().length != upperBounds.value().length || upperBounds.value().length != lowerBounds.value().length) {
           throw new IllegalArgumentException("The calibrations, upper bounds, and lower bounds must have the same length");
        }

        this.calibrationTaxa = calibrationTaxa;
        this.upperBounds = upperBounds;
        this.lowerBounds = lowerBounds;
        this.coverage = coverage;
        this.rootFlag = rootFlag;
    }

    @GeneratorInfo(name = "ConditionedMRCAPrior", examples = {"conditionedMRCAPrior.lphy"},
        description = "generate an array of calibrations with optimised parameters for the distributions on the MRCA node.")
    @Override
    public RandomVariable<Calibration[]> sample() {
        // get the parameters
        String[][] calibrationTaxa = getCalibrationTaxa().value();
        Number[] upperBounds = getUpperBounds().value();
        Number[] lowerBounds = getLowerBounds().value();
        double coverage = 0.9;
        if (getCoverage() != null) {
            coverage = getCoverage().value().doubleValue();
        }
        boolean rootFlag = false;
        if (getRootFlag() != null) {
            rootFlag = getRootFlag().value();
        }
        int n = calibrationTaxa.length;

        /*
            step 1: get maximal clade calibrations, figure out who's whose parent
         */
        // build clade calibrations first
        List<Calibration> calibrations = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            Calibration cali = new Calibration(calibrationTaxa[i]);
            calibrations.add(cali);
        }

        // convert a simple array that telling calibrations relationship
        int[] parent = computeParents(rootFlag, calibrations);

        /*
            step 2 : classify nodes, see if it is beta nodes
            beta node is nested within a larger clade along the main branch
            they're using beta*lognormal
            if there is overlapping of bound between the child calibration and the parent then set true
            otherwise it's still not beta node
         */
        boolean[] is_beta_node = mapBetaNodes(n, rootFlag, parent, upperBounds, lowerBounds);

        /*
            step 3: compute the log-targets
            use the methods from CalibrationPrior beast class
         */
        double[] mu = new double[n];
        double[] sigma2 = new double[n];
        for (int i = 0; i < n; i++) {
            double lo = lowerBounds[i].doubleValue();
            double hi = upperBounds[i].doubleValue();
            mu[i]     = computeLogTargetsMu(lo, hi, coverage);
            sigma2[i] = computeLogTargetsSigma2(lo, hi, coverage);
        }

        /*
            step 4: compute edges for beta nodes
         */
        // count edges <- is_beta_node and not root
        int nEdges = getEdgesNumber(n, is_beta_node, parent);

        // allocate outputs
        int[] edges = new int[nEdges];
        String[] edgeNames = new String[nEdges];

        double[][] A_mean = new double[nEdges][nEdges]; // default 0.0
        double[] b_mean = new double[nEdges];
        double[] b_var  = new double[nEdges];

        // fill
        int idx = 0;
        for (int i = 0; i < n; i++) {
            if (is_beta_node[i] && parent[i] != -1) {
                int par_i = parent[i];

                edges[idx] = i;
                edgeNames[idx] = par_i + "_" + i;   // parent_child

                double m_edge = mu[i] - mu[par_i];
                double v_edge = sigma2[i];       // approx relative variance

                b_mean[idx] = m_edge;
                b_var[idx]  = v_edge;
                A_mean[idx][idx] = 1.0;

                idx++;
            }
        }

        /*
            step 5: solve for beta parameters
         */

        java.util.Map<String, Double> edgeAlpha = new java.util.HashMap<>();
        java.util.Map<String, Double> edgeBeta  = new java.util.HashMap<>();
        extractBetaParams(nEdges, b_mean, b_var, edgeAlpha, edgeNames, edgeBeta);

        /*
            step 6: simulation
         */
        double[] W = new double[n];
        calculateNodeAges(n, rootFlag, W, mu, sigma2, parent, is_beta_node, edgeAlpha, edgeBeta);

        /*
            step 7: map output
         */
        Calibration[] calibrationOuts = new Calibration[n];
        for (int i = 0; i < n; i++) {
            Calibration calibration = new Calibration(calibrations.get(i).getTaxa(), W[i]);
            calibrationOuts[i] = calibration;

        }

        return new RandomVariable<>("", calibrationOuts, this);
    }

    private static void calculateNodeAges(int n, boolean rootFlag, double[] W, double[] mu, double[] sigma2, int[] parent, boolean[] is_beta_node, Map<String, Double> edgeAlpha, Map<String, Double> edgeBeta) {
        // mark unassigned
        for (int i = 0; i < n; i++) {
            W[i] = Double.NaN;
        }

        boolean[] inStack = new boolean[n];

        for (int i = 0; i < n; i++) {
            calculateByOrder(i, rootFlag, W, mu, sigma2, parent, is_beta_node, edgeAlpha, edgeBeta, inStack);
        }
    }

    private static void calculateByOrder(int i, boolean rootFlag, double[] W, double[] mu, double[] sigma2, int[] parent, boolean[] isBetaNode, Map<String, Double> edgeAlpha, Map<String, Double> edgeBeta, boolean[] inStack) {
        // already assigned
        if (Double.isFinite(W[i])) {
            return;
        }

        // detect cycle
        if (inStack[i]) {
            throw new IllegalStateException("Cycle detected at node " + i);
        }

        inStack[i] = true;

        int p = parent[i];

        NormalDistribution nd = new NormalDistribution(0, 1);

        // root or no parent → plain lognormal
        if (p == -1 || (rootFlag && i == 0)) {
            W[i] = Math.exp(mu[i] + Math.sqrt(sigma2[i]) * nd.sample());
            inStack[i] = false;
            return;
        }

        // ensure parent computed first
        calculateByOrder(p, rootFlag, W, mu, sigma2, parent, isBetaNode, edgeAlpha, edgeBeta, inStack);

        if (!Double.isFinite(W[p]) || W[p] <= 0.0) {
            throw new IllegalStateException("Parent age invalid: parent=" + p);
        }

        if (!isBetaNode[i]) {
            TruncatedLogNormal tln = new TruncatedLogNormal(new Value<>("", mu[i]), new Value<>("", Math.sqrt(sigma2[i])), new Value<>("", W[p]));
            W[i] = tln.sample().value();
        } else {
            String key = p + "_" + i;
            double a = edgeAlpha.get(key);
            double b = edgeBeta.get(key);
            BetaDistribution bd = new BetaDistribution(a, b);
            W[i] = W[p] * bd.sample();
        }

        inStack[i] = false;
    }

    private static void extractBetaParams(int nEdges, double[] b_mean, double[] b_var, Map<String, Double> edgeAlpha, String[] edgeNames, Map<String, Double> edgeBeta) {
        // m_hat <- solve(A_mean, b_mean)
        double[] m_hat = new double[nEdges];
        for (int j = 0; j < nEdges; j++) {
            m_hat[j] = b_mean[j];   // A_mean is identity
        }

        // v_hat <- nnls(A_mean, b_var)$x
        double[] v_hat = new double[nEdges];
        for (int j = 0; j < nEdges; j++) {
            v_hat[j] = b_var[j];    // A_mean is identity
        }

        // v_hat <- pmax(v_hat, 1e-6)
        for (int j = 0; j < nEdges; j++) {
            if (v_hat[j] < 1e-6) {
                v_hat[j] = 1e-6;
            }
        }

        // beta_params <- lapply(...)


        for (int j = 0; j < nEdges; j++) {
            double[] ab = invertLogMomentsToBeta(m_hat[j], v_hat[j]);
            edgeAlpha.put(edgeNames[j], ab[0]);
            edgeBeta.put(edgeNames[j], ab[1]);
        }
    }

    private static int getEdgesNumber(int n, boolean[] is_beta_node, int[] parent) {
        int nEdges = 0;
        for (int i = 0; i < n; i++) {
            if (is_beta_node[i] && parent[i] != -1) {
                nEdges++;
            }
        }
        return nEdges;
    }

    public static boolean[] mapBetaNodes(int n, boolean rootFlag, int[] parent, Number[] upperBounds, Number[] lowerBounds) {
        boolean[] is_beta_node = new boolean[n];
        Arrays.fill(is_beta_node, true);

        for (int i = 0; i < n; i++) {
            if (rootFlag && i==0) {
                is_beta_node[i] = false; // root is always lognormal
                continue;
            }

            int par = parent[i];
            if (par != -1) {
                // if child interval does not overlap parent, treat as new root
                if (upperBounds[i].doubleValue() <= lowerBounds[par].doubleValue()) {
                    is_beta_node[i] = false;
                }
            } else {
                is_beta_node[i] = false;
            }
        }
        return is_beta_node;
    }

    // public for unit test
    public static int[] computeParents(boolean rootFlag, List<Calibration> calibrations) {
        int n = calibrations.size();
        // check if the first calibration is the root, otherwise no root
        final int rootIdx = rootFlag ? 0 : -1;

        int[] parent = new int[n];
        Arrays.fill(parent, -1);

        // get taxa names length
        int[] size = new int[n];
        for (int i = 0; i < n; i++) {
            size[i] = calibrations.get(i).getTaxa().length;
        }

        for (int i = 0; i < n; i++) {

            // if root, then assign its parent to NA
            if (rootFlag && i == rootIdx) {
                parent[i] = -1;
                continue;
            }

            int best = -1;  // best parent candidate index found so far
            int bestSize = Integer.MAX_VALUE;

            // try every calibration j as a potential parent of i.
            for (int j = 0; j < n; j++) {
                if (j == i) continue; //skip itself
                Calibration candidate = calibrations.get(j);

                // check if there are super sets of it
                boolean[] sup = isSuperSetOf(candidate, calibrations);
                if (!sup[i]) continue;
                // parent should be larger than the child clade
                if (size[j] <= size[i]) continue;

                // among all supersets, pick the one with the smallest taxa count (tightest enclosing clade).
                if (size[j] < bestSize) {
                    best = j;
                    bestSize = size[j];
                }
            }

            if (best != -1) {
                parent[i] = best;
                continue;
            }

            // if no parent was found:
            // - has root calibration: attach top-level calibrations to root
            // - no root calibration: keep as NA
            parent[i] = rootFlag ? rootIdx : -1;
        }

        return parent;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(calibrationTaxaName, calibrationTaxa);
        map.put(upperBoundName, upperBounds);
        map.put(lowerBoundName, lowerBounds);
        if (rootFlag != null) map.put(rootFlagName, rootFlag);
        if (coverage != null) map.put(coverageName, coverage);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(calibrationTaxaName)) calibrationTaxa = value;
        else if (paramName.equals(rootFlagName)) rootFlag = value;
        else if (paramName.equals(upperBoundName)) upperBounds = value;
        else if (paramName.equals(lowerBoundName)) lowerBounds = value;
        else if (paramName.equals(coverageName)) coverage = value;
        else throw new IllegalArgumentException("Unknown param: " + paramName);
    }

    public Value<String[][]> getCalibrationTaxa() {
        return calibrationTaxa;
    }
    public Value<Number[]> getUpperBounds() {
        return upperBounds;
    }
    public Value<Number[]> getLowerBounds() {
        return lowerBounds;
    }
    public Value<Number> getCoverage() {
        return coverage;
    }
    public Value<Boolean> getRootFlag() {
        return rootFlag;
    }
}
