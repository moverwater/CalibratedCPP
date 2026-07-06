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

public class ConditionedMRCAPrior implements GenerativeDistribution<CalibrationArray> {

    public static final String calibrationsParamName = "calibrations";
    public static final String coverageName = "p";

    Value<?> calibrationsInput;
    Value<Number> coverage;

    public ConditionedMRCAPrior(
            @ParameterInfo(name = calibrationsParamName, description = "array of calibration constraints, each created with calibration(taxa, upper, lower, root?)") Value<?> calibrationsInput,
            @ParameterInfo(name = coverageName, description = "confidence level for probability mass after truncation, default 0.9", optional = true) Value<Number> coverage) {
        if (calibrationsInput == null)
            throw new IllegalArgumentException("calibrations must be provided");
        Calibration[] specs = extractCalibrations(calibrationsInput);
        if (specs.length < 1)
            throw new IllegalArgumentException("At least one calibration must be provided");
        for (Calibration c : specs) {
            if (c.getUpper() == null || c.getLower() == null)
                throw new IllegalArgumentException("Each calibration must have upper and lower bounds — use calibration(taxa=..., upper=..., lower=...)");
        }
        this.calibrationsInput = calibrationsInput;
        this.coverage = coverage;
    }

    @GeneratorInfo(name = "ConditionedMRCAPrior", examples = {"conditionedMRCAPrior.lphy"},
            description = "Generates an array of calibrated node ages with optimised distribution parameters for the MRCA nodes.")
    @Override
    public RandomVariable<CalibrationArray> sample() {
        Calibration[] calibrationSpecs = extractCalibrations(calibrationsInput);
        int n = calibrationSpecs.length;

        Double[] upperBounds = new Double[n];
        Double[] lowerBounds = new Double[n];
        for (int i = 0; i < n; i++) {
            upperBounds[i] = calibrationSpecs[i].getUpper();
            lowerBounds[i] = calibrationSpecs[i].getLower();
        }

        double p = 0.9;
        if (coverage != null) p = coverage.value().doubleValue();

        List<Calibration> calibrationList = new ArrayList<>(Arrays.asList(calibrationSpecs));
        int[] parent = computeParents(calibrationList);
        boolean[] is_beta_node = mapBetaNodes(calibrationSpecs, parent);

        double[] mu = new double[n];
        double[] sigma2 = new double[n];
        for (int i = 0; i < n; i++) {
            mu[i]     = computeLogTargetsMu(lowerBounds[i], upperBounds[i], p);
            sigma2[i] = computeLogTargetsSigma2(lowerBounds[i], upperBounds[i], p);
        }

        int nEdges = getEdgesNumber(n, is_beta_node, parent);
        int[] edges = new int[nEdges];
        String[] edgeNames = new String[nEdges];
        double[][] A_mean = new double[nEdges][nEdges];
        double[] b_mean = new double[nEdges];
        double[] b_var  = new double[nEdges];

        int idx = 0;
        for (int i = 0; i < n; i++) {
            if (is_beta_node[i] && parent[i] != -1) {
                int par_i = parent[i];
                edges[idx] = i;
                edgeNames[idx] = par_i + "_" + i;
                b_mean[idx] = mu[i] - mu[par_i];
                b_var[idx]  = sigma2[i] - sigma2[par_i];
                A_mean[idx][idx] = 1.0;
                idx++;
            }
        }

        Map<String, Double> edgeAlpha = new HashMap<>();
        Map<String, Double> edgeBeta  = new HashMap<>();
        extractBetaParams(nEdges, b_mean, b_var, edgeAlpha, edgeNames, edgeBeta);

        double[] W = new double[n];
        calculateNodeAges(n, W, mu, sigma2, parent, is_beta_node, edgeAlpha, edgeBeta);

        Calibration[] result = new Calibration[n];
        for (int i = 0; i < n; i++) {
            result[i] = new Calibration(calibrationSpecs[i].getTaxa(), W[i]);
        }
        return new RandomVariable<>("", new CalibrationArray(result), this);
    }

    private static void calculateNodeAges(int n, double[] W, double[] mu, double[] sigma2, int[] parent,
                                          boolean[] is_beta_node, Map<String, Double> edgeAlpha, Map<String, Double> edgeBeta) {
        for (int i = 0; i < n; i++) W[i] = Double.NaN;
        boolean[] inStack = new boolean[n];
        for (int i = 0; i < n; i++) {
            calculateByOrder(i, W, mu, sigma2, parent, is_beta_node, edgeAlpha, edgeBeta, inStack);
        }
    }

    private static void calculateByOrder(int i, double[] W, double[] mu, double[] sigma2, int[] parent,
                                         boolean[] isBetaNode, Map<String, Double> edgeAlpha, Map<String, Double> edgeBeta, boolean[] inStack) {
        if (Double.isFinite(W[i])) return;
        if (inStack[i]) throw new IllegalStateException("Cycle detected at node " + i);
        inStack[i] = true;

        int p = parent[i];
        NormalDistribution nd = new NormalDistribution(0, 1);

        if (p == -1) {
            W[i] = Math.exp(mu[i] + Math.sqrt(sigma2[i]) * nd.sample());
            inStack[i] = false;
            return;
        }

        calculateByOrder(p, W, mu, sigma2, parent, isBetaNode, edgeAlpha, edgeBeta, inStack);

        if (!Double.isFinite(W[p]) || W[p] <= 1e-8) {
            W[i] = (0.5 + Math.random() * 0.4) * Math.max(W[p], 1e-6);
            inStack[i] = false;
            return;
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

    private static void extractBetaParams(int nEdges, double[] b_mean, double[] b_var,
                                          Map<String, Double> edgeAlpha, String[] edgeNames, Map<String, Double> edgeBeta) {
        double[] m_hat = Arrays.copyOf(b_mean, nEdges);
        double[] v_hat = Arrays.copyOf(b_var, nEdges);
        for (int j = 0; j < nEdges; j++) {
            if (v_hat[j] <= 0) v_hat[j] = 1e-8;
        }
        for (int j = 0; j < nEdges; j++) {
            double[] ab = invertLogMomentsToBeta(m_hat[j], v_hat[j]);
            edgeAlpha.put(edgeNames[j], ab[0]);
            edgeBeta.put(edgeNames[j], ab[1]);
        }
    }

    private static int getEdgesNumber(int n, boolean[] is_beta_node, int[] parent) {
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (is_beta_node[i] && parent[i] != -1) count++;
        }
        return count;
    }

    // public for unit tests
    public static boolean[] mapBetaNodes(Calibration[] calibrations, int[] parent) {
        int n = calibrations.length;
        boolean[] is_beta_node = new boolean[n];
        Arrays.fill(is_beta_node, true);

        for (int i = 0; i < n; i++) {
            int par = parent[i];
            if (par != -1) {
                if (calibrations[i].getUpper() <= calibrations[par].getLower()) {
                    is_beta_node[i] = false;
                }
            } else {
                is_beta_node[i] = false;
            }
        }
        return is_beta_node;
    }

    // public for unit tests
    public static int[] computeParents(List<Calibration> calibrations) {
        int n = calibrations.size();

        int[] parent = new int[n];
        Arrays.fill(parent, -1);

        int[] size = new int[n];
        for (int i = 0; i < n; i++) {
            size[i] = calibrations.get(i).getTaxa().length;
        }

        for (int i = 0; i < n; i++) {
            int best = -1;
            int bestSize = Integer.MAX_VALUE;

            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                Calibration candidate = calibrations.get(j);
                boolean[] sup = isSuperSetOf(candidate, calibrations);
                if (!sup[i]) continue;
                if (size[j] <= size[i]) continue;
                if (size[j] < bestSize) {
                    best = j;
                    bestSize = size[j];
                }
            }

            parent[i] = best; // -1 if no superset found (root of its tree)
        }
        return parent;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(calibrationsParamName, calibrationsInput);
        if (coverage != null) map.put(coverageName, coverage);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(calibrationsParamName)) calibrationsInput = value;
        else if (paramName.equals(coverageName)) coverage = value;
        else throw new IllegalArgumentException("Unknown param: " + paramName);
    }

    /** Returns the input calibration specs (with taxa, upper, lower) as a properly-typed array. */
    public Value<Calibration[]> getCalibrations() {
        return new Value<>("", extractCalibrations(calibrationsInput));
    }

    public Value<Number> getCoverage() {
        return coverage;
    }
}
