package calibratedcpp.lphy.prior;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;

public class ConditionedPriorUtils {

    // methods copied from calibratedcpp-beast/calibrationprior/CalibrationPrior.java
    public static double computeLogTargetsMu(double tLo, double tHi, double p){
        // Replaced NormalDistribution.inverseCumulativeProbability with Erf.erfInv
        // This calculates the z-score for the given coverage probability
        double z = Math.sqrt(2.0) * Erf.erfInv(p);
        double sigma = (Math.log(tHi) - Math.log(tLo)) / (2 * z);
        return Math.log(tLo) + z * sigma;
    }

    public static double computeLogTargetsSigma2(double tLo, double tHi, double p){
        // Replaced NormalDistribution.inverseCumulativeProbability with Erf.erfInv
        // This calculates the z-score for the given coverage probability
        double z = Math.sqrt(2.0) * Erf.erfInv(p);
        double sigma = (Math.log(tHi) - Math.log(tLo)) / (2 * z);
        return sigma * sigma;
    }

    public static double[] invertLogMomentsToBeta(double m, double v) {
        double Ey = Math.exp(m + v / 2);
        double Vy = Math.exp(2 * m + v) * (Math.exp(v) - 1);
        if (Vy <= 0) Vy = 1e-8;
        Ey = Math.min(1 - 1e-8, Math.max(1e-8, Ey));

        double common = Ey * (1 - Ey) / Vy - 1;
        double a = Math.max(1e-3, Ey * common);
        double b = Math.max(1e-3, (1 - Ey) * common);

        for (int iter = 0; iter < 40; iter++) {
            double psiA = Gamma.digamma(a);
            double psiAB = Gamma.digamma(a + b);
            double f1 = psiA - psiAB - m;
            double f2 = Gamma.trigamma(a) - Gamma.trigamma(a + b) - v;
            if (Math.abs(f1) < 1e-8 && Math.abs(f2) < 1e-8) break;
            a -= 0.5 * f1;
            b -= 0.5 * f2;
            a = Math.max(a, 1e-6);
            b = Math.max(b, 1e-6);
        }
        return new double[]{a, b};
    }
}
