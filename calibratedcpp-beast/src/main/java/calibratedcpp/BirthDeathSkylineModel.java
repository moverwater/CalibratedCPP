package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;

import java.util.*;

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
                specifiedRates++;
                if (!whichSpecified.isEmpty()) whichSpecified.append(", ");
                whichSpecified.append(rateInputs[i].getName());

                // --- CHANGED LOGIC HERE ---
                // Only enforce dimension check IF change times are explicitly provided.
                // If timeP is null, we allow any dimension for rateP (implies equidistant).
                if (timeP != null) {
                    int rateDim = rateP.getDimension();
                    int timeDim = timeP.getDimension();
                    if (rateDim != timeDim + 1) {
                        throw new IllegalArgumentException("Dimension mismatch for " + rateInputs[i].getName() +
                                ": Explicit change times provided ("+timeDim+"), so rate dimension must be ("+(timeDim+1)+"). Found: " + rateDim);
                    }
                }
            }
        }

        if (reproductiveNumberInput.get() != null && turnoverInput.get() != null)
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover.");

        if (rhoInput.get() == null)
            throw new IllegalArgumentException("Sampling probability (rho) must be specified.");

        if (specifiedRates != 2)
            throw new IllegalArgumentException("Exactly TWO rates must be specified. Found " + specifiedRates + " (" + whichSpecified + ")");

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
        double rootHeight = treeInput.get().getRoot().getHeight();
        double originVal = !conditionOnRoot ? originInput.get().getValue() : 0.0;
        double maxT = conditionOnRoot ? rootHeight : originVal;

        if (!conditionOnRoot && rootHeight >= originVal) return false;

        // 1. Process Times for each parameter individually
        // This handles the "Equidistant" generation if inputs are null
        // --- Internal State ---
        // We store the processed change times for each parameter individually
        List<Double> bTimes = processInput(birthRateInput, birthRateChangeTimesInput, relativeFlags[0], reverseFlags[0], maxT);
        List<Double> dTimes = processInput(deathRateInput, deathRateChangeTimesInput, relativeFlags[1], reverseFlags[1], maxT);
        List<Double> rTimes = processInput(reproductiveNumberInput, reproductiveNumberChangeTimesInput, relativeFlags[2], reverseFlags[2], maxT);
        List<Double> divTimes = processInput(diversificationRateInput, diversificationRateChangeTimesInput, relativeFlags[3], reverseFlags[3], maxT);
        List<Double> tTimes = processInput(turnoverInput, turnoverChangeTimesInput, relativeFlags[4], reverseFlags[4], maxT);

        // 2. Create Master Timeline (Union of all times)
        SortedSet<Double> timesSet = new TreeSet<>();
        timesSet.add(0.0);
        timesSet.addAll(bTimes); timesSet.addAll(dTimes);
        timesSet.addAll(rTimes); timesSet.addAll(divTimes);
        timesSet.addAll(tTimes);

        intervalStartTimes = timesSet.stream().mapToDouble(d -> d).toArray();
        int n = intervalStartTimes.length;
        lambda = new double[n];
        r = new double[n];
        cumulativeIntegral = new double[n];
        cumulativeExpR = new double[n];

        double logRunningSum = Double.NEGATIVE_INFINITY;
        double rRunningSum = 0.0;
        cumulativeIntegral[0] = Double.NEGATIVE_INFINITY;
        cumulativeExpR[0] = 0.0;

        // 3. Iterate through Master Intervals
        for (int j = 0; j < n; j++) {
            double t = intervalStartTimes[j];

            // Lookup the rate for this specific time 't' using that parameter's specific time list
            double vL = getVal(birthRateInput.get(), bTimes, t, reverseFlags[0]);
            double vM = getVal(deathRateInput.get(), dTimes, t, reverseFlags[1]);
            double vR = getVal(reproductiveNumberInput.get(), rTimes, t, reverseFlags[2]);
            double vD = getVal(diversificationRateInput.get(), divTimes, t, reverseFlags[3]);
            double vT = getVal(turnoverInput.get(), tTimes, t, reverseFlags[4]);

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

    private List<Double> processInput(Input<RealParameter> rateInput, Input<RealParameter> timeInput,
                                      boolean relative, boolean reverse, double maxTime) {
        List<Double> times = new ArrayList<>();
        if (rateInput.get() == null) return times;

        RealParameter rateP = rateInput.get();
        RealParameter timeP = timeInput.get();

        if (timeP != null) {
            // --- Explicit Times ---
            Double[] vals = timeP.getValues();
            Arrays.sort(vals); // Ensure sorted
            for (double v : vals) {
                double tAbs = relative ? (reverse ? maxTime * (1.0 - v) : maxTime * v)
                        : (reverse ? maxTime - v : v);
                if (tAbs > 1e-10 && tAbs < maxTime - 1e-10) times.add(tAbs);
            }
        } else {
            // --- Implicit Equidistant Times ---
            // If rate has dimension K, we need K-1 cut points.
            int numIntervals = rateP.getDimension();
            if (numIntervals > 1) {
                double width = maxTime / numIntervals;
                for (int i = 1; i < numIntervals; i++) {
                    // Equidistant splits: 1/K, 2/K, etc.
                    // We don't care about reverse here because symmetric equidistant splits are the same reversed.
                    times.add(width * i);
                }
            }
        }
        Collections.sort(times);
        return times;
    }

    private double getVal(RealParameter p, List<Double> cuts, double t, boolean rev) {
        if (p == null) return 0;

        // Binary search to find which interval 't' falls into
        int idx = Collections.binarySearch(cuts, t);

        // binarySearch returns:
        // >= 0: Exact match found. We move to the NEXT interval (right-continuous).
        // < 0:  (-insertionPoint - 1). Insertion point is the index of the first element greater than key.

        int intervalIndex;
        if (idx >= 0) {
            intervalIndex = idx + 1;
        } else {
            intervalIndex = -idx - 1;
        }

        // Map interval index to parameter array index
        // If Reverse (Present -> Past): Rate[0] is most recent (time 0).
        // If Forward (Past -> Present): Rate[Dim-1] is most recent.

        int pIdx = rev ? intervalIndex : (p.getDimension() - 1 - intervalIndex);

        // Safety clamp
        return p.getValue(Math.max(0, Math.min(pIdx, p.getDimension() - 1)));
    }

    // --- Math Helpers (unchanged) ---

    private double calculateSegment(double l, double r_val, double dt, double expOffset) {
        double logTerm;
        double x = r_val * dt;
        if (Math.abs(r_val) < 1e-9) {
            logTerm = Math.log(rho) + Math.log(l) + Math.log(dt);
        } else if (r_val > 0) {
            double logDiff = Math.log(Math.expm1(x));
            logTerm = Math.log(rho) + Math.log(l) - Math.log(r_val) + logDiff;
        } else {
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

    private double logSumExp(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b; if (b == Double.NEGATIVE_INFINITY) return a;
        return Math.max(a, b) + Math.log1p(Math.exp(-Math.abs(a - b)));
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