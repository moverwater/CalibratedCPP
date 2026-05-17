package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Gamma;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;

/**
 * @author Marcus Overwater
 */
@Description("Implementation of the Calibrated Coalescent Point Process where individual lifetimes follow an " +
        "Erlang distribution (integer-shape Gamma) and births happen at a constant rate.")
public class CalibratedAgeDependentBirthDeathModel extends CalibratedCoalescentPointProcess {

    Input<ParametricDistribution> lifetimeDistributionInput = new Input<>("lifetimeDistribution",
            "Distribution of the lifetime of an individual.");
    Input<RealParameter> birthRateInput = new Input<>("birthRate", "The birth rate.");

    protected boolean lifetimesAreErlang;
    protected double birthRate;
    protected double rho;
    protected Gamma gammaDistribution;
    protected int n;        // Erlang shape (positive integer)
    protected double theta; // Erlang rate (1/scale)

    // Precomputed for F(t) = gammaConst + sum_j alphas[j] * exp(roots[j] * t)
    protected Complex[] roots;
    protected Complex[] alphas;
    protected double gammaConst;

    @Override
    public void initAndValidate() {
        lifetimesAreErlang = lifetimeDistributionInput.get() instanceof Gamma;
        super.initAndValidate();
    }

    @Override
    public void updateModel() {
        super.updateModel();
        preCalc();
    }

    public void preCalc() {
        birthRate = birthRateInput.get().getValue();

        if (!lifetimesAreErlang) return;

        gammaDistribution = (Gamma) lifetimeDistributionInput.get();
        double shapeParam = gammaDistribution.alphaInput.get().getArrayValue();
        n = (int) Math.round(shapeParam);
        if (Math.abs(shapeParam - n) > 1e-10 || n < 1) {
            throw new IllegalArgumentException(
                    "Shape parameter must be a positive integer for Erlang distribution, got " + shapeParam);
        }
        theta = 1.0 / gammaDistribution.betaInput.get().getArrayValue();

        // R_n(x) = (x+theta)^n - lambda * [(x+theta)^n - theta^n] / x
        double[] coeffs = buildRnCoefficients(n, theta, birthRate);

        roots = new LaguerreSolver(1e-12).solveAllComplex(coeffs, 1.0);

        // gammaConst = theta^n / R_n(0) = theta^n / coeffs[0]
        gammaConst = Math.pow(theta, n) / coeffs[0];

        // alpha_j = (x_j + theta)^n / (x_j * R_n'(x_j))
        alphas = computeResidues(roots, n, theta, coeffs);
    }

    /**
     * Coefficients of R_n(x) in ascending power order (c[0] = constant, c[n] = leading).
     * Expanding R_n(x) = (x+theta)^n - lambda*[(x+theta)^n - theta^n]/x and collecting powers:
     * c[n] = 1; for m < n: c[m] = C(n,n-m)*theta^(n-m) - lambda*C(n,n-m-1)*theta^(n-m-1).
     * For n=2 this gives Q(x) from Lambert & Stadler (2013) Proposition 6.
     */
    private double[] buildRnCoefficients(int n, double theta, double lambda) {
        double[] c = new double[n + 1];
        c[n] = 1.0;
        for (int m = 0; m < n; m++) {
            int k = n - m;
            c[m] = binomial(n, k) * Math.pow(theta, k)
                 - lambda * binomial(n, k - 1) * Math.pow(theta, k - 1);
        }
        return c;
    }

    private Complex[] computeResidues(Complex[] roots, int n, double theta, double[] coeffs) {
        Complex[] result = new Complex[n];
        for (int j = 0; j < n; j++) {
            Complex xj = roots[j];
            Complex numerator = xj.add(theta).pow(n);
            Complex rPrime = evaluateDerivative(coeffs, xj);
            result[j] = numerator.divide(xj.multiply(rPrime));
        }
        return result;
    }

    private Complex evaluateDerivative(double[] coeffs, Complex x) {
        // sum_{m=1}^{n} m * coeffs[m] * x^{m-1}
        Complex result = Complex.ZERO;
        Complex xPow = Complex.ONE;
        for (int m = 1; m < coeffs.length; m++) {
            result = result.add(xPow.multiply(m * coeffs[m]));
            xPow = xPow.multiply(x);
        }
        return result;
    }

    /**
     * F_p(t) = (1-rho) + rho*gammaConst + rho * Re[sum_j alphas[j] * exp(roots[j]*t)]
     * Complex roots come in conjugate pairs so imaginary parts cancel in the real-valued sum.
     */
    private double computeFp(double t) {
        double sum = (1.0 - rho) + rho * gammaConst;
        for (int j = 0; j < n; j++) {
            sum += rho * alphas[j].multiply(roots[j].multiply(t).exp()).getReal();
        }
        return sum;
    }

    /** F_p'(t) = rho * Re[sum_j alphas[j] * roots[j] * exp(roots[j]*t)] */
    private double computeFpPrime(double t) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += alphas[j].multiply(roots[j]).multiply(roots[j].multiply(t).exp()).getReal();
        }
        return rho * sum;
    }

    @Override
    public double calculateLogNodeAgeDensity(double time) {
        if (!lifetimesAreErlang) return 0.0;
        double fp = computeFp(time);
        double fpPrime = computeFpPrime(time);
        // f_p(t) = F_p'(t) / F_p(t)^2
        return Math.log(fpPrime) - 2.0 * Math.log(fp);
    }

    @Override
    public double calculateLogNodeAgeCDF(double time) {
        if (!lifetimesAreErlang) return 0.0;
        double fp = computeFp(time);
        // P(H_p < t) = 1 - 1/F_p(t) = (F_p(t) - 1) / F_p(t)
        return Math.log(fp - 1.0) - Math.log(fp);
    }

    private static long binomial(int n, int k) {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        k = Math.min(k, n - k);
        long result = 1;
        for (int i = 0; i < k; i++) {
            result = result * (n - i) / (i + 1);
        }
        return result;
    }
}
