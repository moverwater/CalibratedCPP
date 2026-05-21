package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.Calibration;

import java.util.*;

    /*
       The mathematical and sampling methods for CPP
    */
public class CPPUtils {
    // resolve bd rates from diversification and turnover rates
    public static double[] getBD(double d, double t) {
        double lambda = d / (1.0 - t);
        double mu = d * t / (1.0 - t);
        return new double[]{lambda, mu};
    }

    // ****** mathematical methods ******
    public static double CDF(double b, double d, double rho, double t) {
        double p;
        double r = b - d;
        double A = rho * b;
        double B = b * (1- rho) - d;
        if (Math.abs(r) < 1e-10) {
            // when diversification rate ~ 0
            p = A * t / (1 + A * t);
        } else if (r < 0){
            double exp_rt = Math.exp(r * t);
            p = A * (1 - exp_rt) / (-A * exp_rt - B);
        } else{
            // r > 0
            double exp_neg_rt = Math.exp(-r * t);
            p = A * (1-exp_neg_rt) / (A + B * exp_neg_rt);
        }
        return p;
    }

    public static double inverseCDF(double b, double d, double rho, double p) {
        double r = b - d;
        double A = rho * b;
        double B = b * (1 - rho) - d;
        double t;

        if (Math.abs(r) < 1e-10) {
            t = p / (A * (1.0 - p));
        } else if (r < 0) {
            // Subcritical: mirrors CalibratedBirthDeathModel.calculateLogNodeAgeCDF which uses
            // Math.log1p(-exp_rt) and Math.log(-A*exp_rt - B) for numerical stability.
            // A + p*B = |B|*(p_max - p), so log(u) = log|B| + log(p_max-p) - log(A) - log1p(-p)
            // where u = exp(r*t), p_max = A/|B|.
            double absB = -B;             // d - b*(1-rho) > 0 for r < 0
            double pmax = A / absB;
            double logU = Math.log(absB) + Math.log(pmax - p) - Math.log(A) - Math.log1p(-p);
            t = logU / r;
        } else {
            t = Math.log(1 + ((b - d) * p) / (b * rho * (1 - p))) / (b - d);
        }

        return t;
    }

    public static double densityBD(double b, double d, double rho, double time) {
        double density;
        double r = b - d;
        double A = rho * b;
        double B = b * (1- rho) - d;
        if (Math.abs(r) < 1e-10) {
             density = A / ((1.0 + A * time) * (1.0 + A * time));
        } else if (r < 0){
            double exp_rt = Math.exp(r * time);
            density = A* r  * r * r * time / (A * exp_rt + B) * (A * exp_rt + B);
        } else {
            double exp_neg_rt = Math.exp(-r*time);
            density = A * r / Math.pow(A + B * exp_neg_rt, 2) * time;
        }
        return density;
    }

    public static double Qdist(double b, double d, double t, int nSims){
        double p;
        double r = b - d;
        if (Math.abs(r) < 1e-10) {
            p = (b * t) / (1 + b * t);
        } else if (r > 0) {
            p = b * (1 - Math.exp(-r * t)) / (b - d * Math.exp(-r * t));
        } else {
            p = b * (Math.exp(r * t) - 1) / (b* Math.exp(r * t) -d );
        }
        return Math.pow(p, nSims);
    }

    public static double transform(double p, double b, double d, int nSims) {
        double t;
        double r = b - d;
        double u = Math.pow(p, 1.0/nSims);
        if (Math.abs(r) < 1e-10) {
            t = u / (b * (1 - u));
        } else {
           t = Math.log((d * Math.pow(p, (double) 1 / nSims) - b) / (b * (Math.pow(p, (double) 1 / nSims) - 1))) / r;
        }
        return t;
    }


    // ****** time sampling methods ******
    // time sampling methods (with condition time optional and lowerTail optional)
    public static double[] sampleTimes(double birthRate, double deathRate, double samplingProbability, double conditionTime, boolean lowerTail, int nSims) {
        // Calculate the CDF value (Q)
        double Q = CDF(birthRate, deathRate, samplingProbability, conditionTime);

        // Array to store the result
        double[] results = new double[nSims];

        // Generate the samples based on the lowerTail flag
        for (int i = 0; i < nSims; i++) {
            double p;
            if (lowerTail) {
                p = Math.random()*Q;
            } else {
                p = Math.random()*(1-Q) + Q;
            }
            results[i] = inverseCDF(birthRate, deathRate, samplingProbability,p);
        }

        return results;
    }

    public static double[] sampleTimes(double birthRate, double deathRate, double samplingProbability, double conditionTime, int nSims) {
        // Calculate the CDF value (Q)
        double Q = CDF(birthRate, deathRate, samplingProbability, conditionTime);

        // Array to store the result
        double[] results = new double[nSims];

        // Generate the samples based on the lowerTail flag
        for (int i = 0; i < nSims; i++) {
            double p;
            // default sample from [0,Q], lowerTail=True
            p = Math.random()*Q;

            results[i] = inverseCDF(birthRate, deathRate, samplingProbability, p);
        }

        return results;
    }

    public static double[] sampleTimes(double birthRate, double deathRate, double samplingProbability, int nSims) {
        // Calculate the CDF value (Q)
        return sampleTimes(birthRate, deathRate, samplingProbability, 0, Double.POSITIVE_INFINITY, nSims);
    }

    public static double[] sampleTimes(double birthRate, double deathRate, double samplingProbability, double lowerTime, double upperTime, int nSims) {
        double r = birthRate - deathRate;
        double A = samplingProbability * birthRate;
        double B = birthRate * (1 - samplingProbability) - deathRate;

        // True CDF supremum: p_max = A/(-B) < 1 in the subcritical regime (r < 0).
        double pmax = (r < -1e-10) ? A / (-B) : 1.0;

        double Qlower = CDF(birthRate, deathRate, samplingProbability, lowerTime);
        double Qupper = Double.isInfinite(upperTime) ? pmax : CDF(birthRate, deathRate, samplingProbability, upperTime);
        double hiEff  = Double.isInfinite(upperTime) ? lowerTime * 10 + 100 : upperTime;

        double[] times = new double[nSims];

        // Degenerate case: CDF is flat on [lowerTime, upperTime] (both endpoints past the
        // saturation time). No CDF-based inversion can return a t in [lowerTime, upperTime]
        // because every p < p_max maps to t near the saturation boundary, not the requested
        // window. Instead, sample directly from the conditional density on [lowerTime, upperTime].
        // From CalibratedBirthDeathModel.calculateLogNodeAgeDensity for r < 0:
        //   log f(t) ≈ const + r·t  (as exp_rt → 0)
        // so f(t) ∝ exp(r·t) = exp(-|r|·t) — a truncated exponential with rate |r|.
        // Inverting its CDF: t = lowerTime - log1p(-u·scale)/|r|,
        //   scale = -expm1(-|r|·(hi-lo)) = 1 - exp(-|r|·(hi-lo)), ≈ 1 for large |r|·(hi-lo).
        if (r < -1e-10 && Qupper - Qlower < 1e-9) {
            double absR  = -r;
            double scale = -Math.expm1(-absR * (hiEff - lowerTime));
            for (int i = 0; i < nSims; i++) {
                times[i] = lowerTime - Math.log1p(-Math.random() * scale) / absR;
            }
            return times;
        }

        // Normal case: clamp p strictly below p_max so the stable inverseCDF
        // (which uses log(p_max - p), mirroring CalibratedBirthDeathModel's
        // Math.log1p(-exp_rt) and Math.log(-A·exp_rt - B)) always gets a valid argument.
        double eps = 1e-12;
        Qlower = Math.max(eps, Math.min(pmax - 2 * eps, Qlower));
        Qupper = Math.max(Qlower + eps, Math.min(pmax - eps, Qupper));

        for (int i = 0; i < nSims; i++) {
            double p = Math.random() * (Qupper - Qlower) + Qlower;
            times[i] = inverseCDF(birthRate, deathRate, samplingProbability, p);
        }

        return times;
    }

    public static int sampleIndex(double[] weights) {
        // normalize weights
        double sum = 0;
        for (double w : weights) sum += w;
        double[] cdf = new double[weights.length];
        cdf[0] = weights[0] / sum;
        for (int i = 1; i < weights.length; i++) {
            cdf[i] = cdf[i - 1] + weights[i] / sum;
        }

        // generate random number
        double num = Math.random();

        // find index
        for (int i = 0; i < cdf.length; i++) {
            if (num <= cdf[i]) return i;
        }
        return weights.length - 1; // fallback
    }



        // ****** check methods ******
    public static int indexOfMin(List<Double> t) {
        int minIndex = 0;
        double minValue = t.get(0);
        for (int i = 1; i < t.size(); i++) {
            if (t.get(i) < minValue) {
                minValue = t.get(i);
                minIndex = i;
            }
        }
        return minIndex;
    }

    public static List<Integer> checkTrues(boolean[] results) {
        int i = 0;
        List<Integer> indices = new ArrayList<>();
        for (boolean val : results) {
            if (val) {
                indices.add(i);
            }
            i ++;
        }
        return indices;
    }

    // returns list of booleans if clade contains cladeCalibrations.get(i) under the partial order of set inclusion
    public static boolean[] isSuperSetOf(Calibration clade, List<Calibration> cladeCalibrations) {

        Set<String> cladeTaxa = new HashSet<>(Arrays.asList(clade.getTaxa()));

        boolean[] result = new boolean[cladeCalibrations.size()];

        for (int i = 0; i < cladeCalibrations.size(); i++) {

            Calibration entry = cladeCalibrations.get(i);
            Set<String> entryTaxa = new HashSet<>(Arrays.asList(entry.getTaxa()));

            boolean isStrictSubset =
                    entryTaxa.size() < cladeTaxa.size() &&
                            cladeTaxa.containsAll(entryTaxa);

            result[i] = isStrictSubset;

            if (isStrictSubset && clade.getAge() != null) {
                double supersetAge = clade.getAge();
                double subsetAge = entry.getAge();

                if (supersetAge < subsetAge) {
                    throw new IllegalArgumentException(
                            "Superset clade " + Arrays.toString(clade.getTaxa()) +
                                    " has age " + supersetAge +
                                    " which is younger than its subset calibration clade " +
                                    Arrays.toString(entry.getTaxa()) +
                                    " with age " + subsetAge +
                                    ". Please double check the clade ages."
                    );
                }
            }
        }

        return result;
    }

    // returns TRUE iff clade is a STRICT subset of cladeCalibrations.get(i)
    public static boolean[] isSubsetOf(Calibration clade, List<Calibration> cladeCalibrations) {

        Set<String> cladeSet = new HashSet<>(Arrays.asList(clade.getTaxa()));

        boolean[] result = new boolean[cladeCalibrations.size()];

        for (int i = 0; i < cladeCalibrations.size(); i++) {

            Calibration entry = cladeCalibrations.get(i);
            Set<String> entrySet = new HashSet<>(Arrays.asList(entry.getTaxa()));

            // strict subset: clade ⊂ entry
            boolean isStrictSubset =
                    cladeSet.size() < entrySet.size() &&
                            entrySet.containsAll(cladeSet);

            result[i] = isStrictSubset;

            // age check: subset (clade) must not be older than its superset (entry)
            if (isStrictSubset && clade.getAge() != null) {
                double supersetAge = entry.getAge();
                double subsetAge = clade.getAge();

                if (subsetAge > supersetAge) {
                    throw new IllegalArgumentException(
                            "Clade " + Arrays.toString(clade.getTaxa()) +
                                    " has age " + subsetAge +
                                    " which is older than its superset calibration clade " +
                                    Arrays.toString(entry.getTaxa()) +
                                    " with age " + supersetAge +
                                    ". Please double check the clade ages."
                    );
                }
            }
        }

        return result;
    }

    // ****** tree methods ******

    public static double simRandomStem(double birthRate, double deathRate, double greaterThan, int nTaxa){
        double Q = Qdist(birthRate, deathRate, greaterThan, nTaxa);

        double p = Math.random()*(1-Q) + Q;

        double t = transform(p, birthRate, deathRate, nTaxa);
        return t;
    }

    // ****** clade methods ******

    public static List<Calibration> getNestedClades(Calibration clade, List<Calibration> cladeCalibrations) {
        boolean[] isNested = isSuperSetOf(clade,cladeCalibrations);
        List<Integer> indices = checkTrues(isNested);
        List<Calibration> subClades = new ArrayList<>();
        int pointer = 0;

        for (Calibration entry : cladeCalibrations){
            if (indices.contains(pointer)) {
                subClades.add(entry);
            }
            pointer ++;
        }
        return subClades;
    }

    public static List<Calibration> getMaximalCalibrations(List<Calibration> cladeCalibrations) {
        List<Calibration> maximalCalibrations = new ArrayList<>();

        for (int i = 0; i < cladeCalibrations.size(); i++) {
            Calibration current = cladeCalibrations.get(i);
            // check if there's a subset of the calibrations
            boolean[] results = isSubsetOf(current, cladeCalibrations);
            if (checkTrues(results).isEmpty()) {
                maximalCalibrations.add(current);
            }
        }

        return maximalCalibrations;
    }

    public static String makeUnique(String base, Set<String> used) {
        String name = base;
        int k = 2;
        while (used.contains(name)) {
            name = base + "_" + k;
            k++;
        }
        used.add(name);
        return name;
    }


}
