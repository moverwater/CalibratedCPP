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
        if (Vy >= Ey * (1 - Ey)) Vy = 0.99 * Ey * (1 - Ey);

        double common = Ey * (1 - Ey) / Vy - 1;
        double a = Math.max(1e-3, Ey * common);
        double b = Math.max(1e-3, (1 - Ey) * common);

        for (int iter = 0; iter < 100; iter++) {
            double triA  = Gamma.trigamma(a);
            double triAB = Gamma.trigamma(a + b);
            double f1 = Gamma.digamma(a) - Gamma.digamma(a + b) - m;
            double f2 = triA - triAB - v;
            if (Math.abs(f1) < 1e-10 && Math.abs(f2) < 1e-10) break;

            // Analytical Jacobian for f1; numerical (finite-diff) for f2 rows
            // since tetragamma is not in Apache Commons Math.
            double J11 = triA - triAB;          // df1/da
            double J12 = -triAB;                // df1/db
            double ha  = Math.max(a * 1e-5, 1e-9);
            double hb  = Math.max(b * 1e-5, 1e-9);
            double J21 = (Gamma.trigamma(a + ha) - Gamma.trigamma(a + ha + b) - triA + triAB) / ha;
            double J22 = (triAB - Gamma.trigamma(a + b + hb)) / hb;

            double det = J11 * J22 - J12 * J21;
            if (Math.abs(det) < 1e-15) break;

            double da = -(f1 * J22 - f2 * J12) / det;
            double db = -(J11 * f2 - J21 * f1) / det;

            // Backtracking: halve step until both parameters stay positive
            double step = 1.0;
            for (int ls = 0; ls < 20 && (a + step * da < 1e-9 || b + step * db < 1e-9); ls++)
                step *= 0.5;

            a = Math.max(a + step * da, 1e-9);
            b = Math.max(b + step * db, 1e-9);
        }
        return new double[]{a, b};
    }
}
