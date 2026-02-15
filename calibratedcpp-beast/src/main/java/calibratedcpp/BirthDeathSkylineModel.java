package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;
import java.util.TreeSet;

/**
 * @author Marcus Overwater
 */

@Description("Extension of CalibratedCoalescentPointProcess, implements node age density and CDF of " +
        "for the BDSKY model with piecewise constant birth and death rates with incomplete extant sampling.")
public class BirthDeathSkylineModel extends CalibratedCoalescentPointProcess {
    public Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "lambda", (RealParameter) null);
    public Input<RealParameter> deathRateInput =
            new Input<>("deathRate", "mu", (RealParameter) null);
    public Input<RealParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "R = lambda / mu", (RealParameter) null);
    public Input<RealParameter> diversificationRateInput =
            new Input<>("diversificationRate", "r = lambda - mu", (RealParameter) null);
    public Input<RealParameter> turnoverInput =
            new Input<>("turnover", "epsilon = mu / lambda", (RealParameter) null);
    public Input<RealParameter> rhoInput =
            new Input<>("rho", "Sampling probability (rho)", (RealParameter) null);

    public Input<RealParameter> birthRateChangeTimesInput =
            new Input<>("birthRateChangeTimes", "Change times for lambda", (RealParameter) null);
    public Input<RealParameter> deathRateChangeTimesInput =
            new Input<>("deathRateChangeTimes", "Change times for mu", (RealParameter) null);
    public Input<RealParameter> reproductiveNumberChangeTimesInput =
            new Input<>("reproductiveNumberChangeTimes", "Change times for R", (RealParameter) null);
    public Input<RealParameter> diversificationRateChangeTimesInput =
            new Input<>("diversificationRateChangeTimes", "Change times for r", (RealParameter) null);
    public Input<RealParameter> turnoverChangeTimesInput =
            new Input<>("turnoverChangeTimes", "Change times for turnover", (RealParameter) null);

    public Input<Boolean> birthRateTimesRelativeInput =
            new Input<>("birthRateTimesRelative", "Relative to height", false);
    public Input<Boolean> deathRateTimesRelativeInput =
            new Input<>("deathRateTimesRelative", "Relative to height", false);
    public Input<Boolean> reproductiveNumberTimesRelativeInput =
            new Input<>("reproductiveNumberTimesRelative", "Relative", false);
    public Input<Boolean> diversificationRateTimesRelativeInput =
            new Input<>("diversificationRateTimesRelative", "Relative", false);
    public Input<Boolean> turnoverTimesRelativeInput =
            new Input<>("turnoverTimesRelative", "Relative", false);

    public Input<BooleanParameter> reverseTimeArraysInput =
            new Input<>("reverseTimeArrays", "True if parameter values are specified from present to root: [birthRate, deathRate, reproductionNumber, diversificationRate, turnover]", (BooleanParameter) null);

    protected double[] intervalStartTimes, lambda, r, cumulativeIntegral, cumulativeExpR;
    protected double rho;
    private boolean[] reverseFlags;
    private boolean[] relativeFlags;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // 1. Group inputs for iterative validation
        Input<RealParameter>[] rateInputs = new Input[]{birthRateInput, deathRateInput, reproductiveNumberInput, diversificationRateInput, turnoverInput};
        Input<RealParameter>[] timeInputs = new Input[]{birthRateChangeTimesInput, deathRateChangeTimesInput, reproductiveNumberChangeTimesInput, diversificationRateChangeTimesInput, turnoverChangeTimesInput};
        Input<Boolean>[] relInputs = new Input[]{birthRateTimesRelativeInput, deathRateTimesRelativeInput, reproductiveNumberTimesRelativeInput, diversificationRateTimesRelativeInput, turnoverTimesRelativeInput};

        relativeFlags = new boolean[5];
        int specifiedRates = 0;
        StringBuilder whichSpecified = new StringBuilder();

        for (int i = 0; i < 5; i++) {
            relativeFlags[i] = relInputs[i].get();
            RealParameter rateP = rateInputs[i].get();
            RealParameter timeP = timeInputs[i].get();

            if (rateP != null) {
                // Check Dimensions
                if (rateInputs[i].get() != null) {
                    int rateDim = rateInputs[i].get().getDimension();
                    int timeDim = (timeInputs[i].get() != null) ? timeInputs[i].get().getDimension() : 0;

                    if (rateDim != timeDim + 1) {
                        throw new IllegalArgumentException("Dimension mismatch for " + rateInputs[i].getName() +
                                ": found " + rateDim + " rates but " + timeDim + " change times. " +
                                "Number of rates must be change times + 1.");
                    }
                }

                // Track specified rates for the "Exactly Two" check
                specifiedRates++;
                if (!whichSpecified.isEmpty()) whichSpecified.append(", ");
                whichSpecified.append(rateInputs[i].getName());

                // Check Relative Constraints
                if (relativeFlags[i]) {
                    if (timeP != null) {
                        for (double val : timeP.getDoubleValues()) {
                            if (val < 0.0 || val > 1.0)
                                throw new IllegalArgumentException("Relative time in " + timeInputs[i].getName() + " must be in [0,1]. Found: " + val);
                        }
                    }
                }
            }
        }

        // 2. Mutual Exclusion & Requirement Checks
        if (reproductiveNumberInput.get() != null && turnoverInput.get() != null)
            throw new IllegalArgumentException("reproductiveNumber and turnover cannot be specified together.");

        if (rhoInput.get() == null)
            throw new IllegalArgumentException("Sampling probability (rho) must be specified.");

        if (specifiedRates != 2)
            throw new IllegalArgumentException("Exactly TWO rates must be specified. Found " + specifiedRates + " (" + whichSpecified + ")");

        // 3. Handle Reverse Flags
        reverseFlags = new boolean[5];
        if (reverseTimeArraysInput.get() != null) {
            BooleanParameter flags = reverseTimeArraysInput.get();
            for (int i = 0; i < Math.min(5, flags.getDimension()); i++) reverseFlags[i] = flags.getValue(i);
        }

        updateIntervals();
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        if (!updateIntervals()) {
            return Double.NEGATIVE_INFINITY;
        }
        return super.calculateTreeLogLikelihood(tree);
    }

    private boolean updateIntervals() {
        rho = rhoInput.get().getValue();

        // Get tree height and origin
        double rootHeight = treeInput.get().getRoot().getHeight();
        double originVal = !conditionOnRoot ? originInput.get().getValue() : 0.0;
        double maxT = conditionOnRoot ? rootHeight : originVal;

        // --- CHECK 1: Root vs Origin ---
        // If we condition on origin, the tree cannot be older than the origin.
        if (!conditionOnRoot && rootHeight >= originVal) {
            return false;
        }

        // --- CHECK 2: Change Times vs Origin ---
        double[] bC = getSorted(birthRateChangeTimesInput.get(), relativeFlags[0], reverseFlags[0], maxT);
        double[] dC = getSorted(deathRateChangeTimesInput.get(), relativeFlags[1], reverseFlags[1], maxT);
        double[] rC = getSorted(reproductiveNumberChangeTimesInput.get(), relativeFlags[2], reverseFlags[2], maxT);
        double[] divC = getSorted(diversificationRateChangeTimesInput.get(), relativeFlags[3], reverseFlags[3], maxT);
        double[] tC = getSorted(turnoverChangeTimesInput.get(), relativeFlags[4], reverseFlags[4], maxT);

        if (isAnyTimeBeyondMax(maxT, bC, dC, rC, divC, tC)) {
            return false;
        }
        // Build the master timeline
        TreeSet<Double> times = new TreeSet<>();
        times.add(0.0);
        for (double t : bC) times.add(t); for (double t : dC) times.add(t);
        for (double t : rC) times.add(t); for (double t : divC) times.add(t);
        for (double t : tC) times.add(t);

        intervalStartTimes = times.stream().mapToDouble(d -> d).toArray();
        int n = intervalStartTimes.length;
        lambda = new double[n]; r = new double[n];
        cumulativeIntegral = new double[n]; cumulativeExpR = new double[n];

        double logRunningSum = Double.NEGATIVE_INFINITY;
        double rRunningSum = 0.0;
        cumulativeIntegral[0] = Double.NEGATIVE_INFINITY;
        cumulativeExpR[0] = 0.0;

        // Rate transformation logic
        for (int j = 0; j < n; j++) {
            double t = intervalStartTimes[j];
            double vL = getVal(birthRateInput.get(), bC, t, reverseFlags[0]);
            double vM = getVal(deathRateInput.get(), dC, t, reverseFlags[1]);
            double vR = getVal(reproductiveNumberInput.get(), rC, t, reverseFlags[2]);
            double vD = getVal(diversificationRateInput.get(), divC, t, reverseFlags[3]);
            double vT = getVal(turnoverInput.get(), tC, t, reverseFlags[4]);

            // Transform user inputs to Lambda and Mu
            double l = 0, m = 0;
            if (birthRateInput.get() != null && deathRateInput.get() != null) { l = vL; m = vM; }
            else if (birthRateInput.get() != null && diversificationRateInput.get() != null) { l = vL; m = vL - vD; }
            else if (birthRateInput.get() != null && reproductiveNumberInput.get() != null) { l = vL; m = vL / vR; }
            else if (birthRateInput.get() != null && turnoverInput.get() != null) { l = vL; m = vL * vT; }
            else if (deathRateInput.get() != null && diversificationRateInput.get() != null) { m = vM; l = vM + vD; }
            else if (deathRateInput.get() != null && turnoverInput.get() != null) { m = vM; l = vM / vT; }
            else if (deathRateInput.get() != null && reproductiveNumberInput.get() != null) { m = vM; l = m * vR;}
            else if (diversificationRateInput.get() != null && reproductiveNumberInput.get() != null) { m = vD / (vR - 1.0); l = m * vR; }
            else if (diversificationRateInput.get() != null && turnoverInput.get() != null) { l = vD / (1.0 - vT); m = l * vT; }

            lambda[j] = l;
            r[j] = l - m;

            if (j < n - 1) {
                double dt = intervalStartTimes[j+1] - t;
                logRunningSum = logSumExp(logRunningSum, calculateSegment(lambda[j], r[j], dt, rRunningSum));
                rRunningSum += (r[j] * dt);
                cumulativeIntegral[j+1] = logRunningSum;
                cumulativeExpR[j+1] = rRunningSum;
            }
        }
        return true;
    }

    private double calculateSegment(double l, double r_val, double dt, double expOffset) {
        double logTerm;
        double x = r_val * dt;

        // --- CASE 1: CRITICAL (r = 0) ---
        if (Math.abs(r_val) < 1e-9) {
            logTerm = Math.log(rho) + Math.log(l) + Math.log(dt);
        }

        // --- CASE 2: SUPER-CRITICAL (r > 0) ---
        else if (r_val > 0) {
            double logDiff = Math.log(Math.expm1(x));
            logTerm = Math.log(rho) + Math.log(l) - Math.log(r_val) + logDiff;
        }

        // --- CASE 3: SUB-CRITICAL (r < 0) ---
        else {
            double logDiff = Math.log(-Math.expm1(x));
            logTerm = Math.log(rho) + Math.log(l) - Math.log(-r_val) + logDiff;
        }

        return logTerm + expOffset;
    }

    @Override
    public double calculateLogNodeAgeDensity(double time) {
        int m = getInterval(time);
        double dtP = time - intervalStartTimes[m];
        double logY = logSumExp(0.0, logSumExp(cumulativeIntegral[m], calculateSegment(lambda[m], r[m], dtP, cumulativeExpR[m])));
        return (Math.log(rho) + Math.log(lambda[m]) + cumulativeExpR[m] + (r[m] * dtP)) - 2.0 * logY;
    }

    @Override
    public double calculateLogNodeAgeCDF(double time) {
        int m = getInterval(time);
        double logInt = logSumExp(cumulativeIntegral[m], calculateSegment(lambda[m], r[m], time - intervalStartTimes[m], cumulativeExpR[m]));
        return logInt - logSumExp(0.0, logInt);
    }

    private int getInterval(double t) {
        int i = Arrays.binarySearch(intervalStartTimes, t);
        return i < 0 ? Math.max(0, -i - 2) : i;
    }

    private boolean isAnyTimeBeyondMax(double maxT, double[]... changePointArrays) {
        for (double[] array : changePointArrays) {
            if (array.length > 0 && array[array.length - 1] > maxT + 1e-10) {
                return true;
            }
        }
        return false;
    }

    private double logSumExp(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b; if (b == Double.NEGATIVE_INFINITY) return a;
        return Math.max(a, b) + Math.log1p(Math.exp(-Math.abs(a - b)));
    }

    private double[] getSorted(RealParameter p, boolean rel, boolean rev, double max) {
        if (p == null || p.getDimension() == 0) return new double[0];
        double[] c = new double[p.getDimension()];
        for (int i = 0; i < c.length; i++) {
            double v = p.getValue(i);
            // Calculate absolute time
            double time = rel ? (rev ? max * (1 - v) : max * v) : (rev ? max - v : v);
            c[i] = time;
        }
        Arrays.sort(c);
        return c;
    }

    private double getVal(RealParameter p, double[] cuts, double t, boolean rev) {
        if (p == null) return 0;

        // Find how many change-points we have passed (starting from t=0)
        int idx = Arrays.binarySearch(cuts, t);
        if (idx < 0) idx = -idx - 1;

        // pIdx logic:
        // If rev = true (Present to Root): idx 0 is the first-rate.
        // If rev = false (Root to Present): We must flip the index.
        int pIdx = rev ? idx : (p.getDimension() - 1 - idx);

        // Safety clamp to prevent index out of bounds
        return p.getValue(Math.max(0, Math.min(pIdx, p.getDimension() - 1)));
    }

    @Override
    public void restore() {
        super.restore();
        updateIntervals();
    }

    @Override
    public boolean requiresRecalculation() {
        super.requiresRecalculation();
        updateIntervals();
        return true;
    }
}
