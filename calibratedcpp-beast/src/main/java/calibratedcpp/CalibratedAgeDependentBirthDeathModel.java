package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Gamma;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;

/**
 * @author Marcus Overwater
 */
@Description("Implementation of the Calibrated Coalescent Point Process where individual lifetimes follow a " +
        "user-specified distribution and births happen at a constant rate. Erlang (integer-shape Gamma) lifetimes " +
        "are handled via closed-form partial fractions; all other distributions use a numerical Volterra IDE solver.")
public class CalibratedAgeDependentBirthDeathModel extends CalibratedCoalescentPointProcess {

    Input<ParametricDistribution> lifetimeDistributionInput = new Input<>("lifetimeDistribution",
            "Distribution of the lifetime of an individual.");
    Input<RealParameter> rhoInput = new Input<>("rho", "Extant sampling probability.");
    Input<RealParameter> birthRateInput = new Input<>("birthRate", "The birth rate.");
    Input<Integer> gridSizeInput = new Input<>("gridSize",
            "Number of grid points for the numerical Volterra IDE solver (used for non-Erlang lifetime distributions).", 1000);

    protected boolean lifetimesAreErlang;
    protected boolean useNumericalSolver;
    protected double birthRate;
    protected double rho;

    // --- Erlang closed-form fields ---
    protected Gamma gammaDistribution;
    protected int n;        // Erlang shape (positive integer)
    protected double theta; // Erlang rate (1/scale)
    protected Complex[] roots;
    protected Complex[] alphas;
    protected double gammaConst;

    // --- Numerical VIDE fields ---
    // gSpline stores G(t) = F(t) - 1, so fp - 1 = rho*G(t) without cancellation
    protected PolynomialSplineFunction gSpline;
    protected PolynomialSplineFunction fPrimeSpline;

    @Override
    public void initAndValidate() {
        ParametricDistribution dist = lifetimeDistributionInput.get();
        lifetimesAreErlang = dist instanceof Gamma;
        useNumericalSolver = (dist != null) && !lifetimesAreErlang;
        super.initAndValidate();
    }

    @Override
    public void updateModel() {
        super.updateModel();
        preCalc();
    }

    public void preCalc() {
        birthRate = birthRateInput.get().getValue();
        rho = rhoInput.get().getValue();

        if (lifetimesAreErlang) {
            preCalcErlang();
        } else if (useNumericalSolver) {
            solveVIDE();
        }
    }

    // -------------------------------------------------------------------------
    // Erlang closed-form path
    // -------------------------------------------------------------------------

    private void preCalcErlang() {
        gammaDistribution = (Gamma) lifetimeDistributionInput.get();
        double shapeParam = gammaDistribution.alphaInput.get().getArrayValue();
        n = (int) Math.round(shapeParam);
        if (Math.abs(shapeParam - n) > 1e-10 || n < 1) {
            throw new IllegalArgumentException(
                    "Shape parameter must be a positive integer for Erlang distribution, got " + shapeParam);
        }
        theta = 1.0 / gammaDistribution.betaInput.get().getArrayValue();

        double[] coeffs = buildRnCoefficients(n, theta, birthRate);
        roots = new LaguerreSolver(1e-12).solveAllComplex(coeffs, 1.0);
        gammaConst = Math.pow(theta, n) / coeffs[0];
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
        Complex result = Complex.ZERO;
        Complex xPow = Complex.ONE;
        for (int m = 1; m < coeffs.length; m++) {
            result = result.add(xPow.multiply(m * coeffs[m]));
            xPow = xPow.multiply(x);
        }
        return result;
    }

    /**
     * Computes scaled exponential sums to prevent overflow for large t.
     *
     * <p>Let {@code maxExp = max_j Re(roots[j]) * t}. Every term is divided by
     * {@code exp(maxExp)} before summing, so the values stay in range regardless of t.
     * The returned array contains:
     * <ul>
     *   <li>[0] maxExp</li>
     *   <li>[1] scaledF  = Re[sum_j alphas[j]          * exp(roots[j]*t - maxExp)]</li>
     *   <li>[2] scaledFP = Re[sum_j alphas[j]*roots[j]  * exp(roots[j]*t - maxExp)]</li>
     * </ul>
     * Then F_p(t) = exp(maxExp) * innerFp  where innerFp = constTerm*exp(-maxExp) + rho*scaledF,
     * and  F_p'(t) = rho * exp(maxExp) * scaledFP.
     * The exp(maxExp) factors cancel in the log density and log CDF formulas.</p>
     */
    private double[] computeScaledSums(double t) {
        double maxExp = 0.0;
        for (int j = 0; j < n; j++) {
            maxExp = Math.max(maxExp, roots[j].getReal() * t);
        }
        double scaledF = 0.0, scaledFP = 0.0;
        for (int j = 0; j < n; j++) {
            Complex expScaled = roots[j].multiply(t).subtract(maxExp).exp();
            scaledF  += alphas[j].multiply(expScaled).getReal();
            scaledFP += alphas[j].multiply(roots[j]).multiply(expScaled).getReal();
        }
        return new double[]{maxExp, scaledF, scaledFP};
    }

    // -------------------------------------------------------------------------
    // Numerical Volterra IDE solver
    // -------------------------------------------------------------------------

    /**
     * Solves F'(t) = lambda*(F(t) - integral_0^t F(t-s)*g(s)ds), F(0) = 1
     * on [0, maxTime] using an implicit trapezoidal scheme.
     *
     * <p>Stores G(t) = F(t) - 1 rather than F(t) directly. Since G(0) = 0 and
     * fp - 1 = rho*G(t), this avoids catastrophic cancellation when evaluating the
     * log CDF at small times where F(t) is close to 1.</p>
     *
     * <p>The VIDE in terms of G(t) is:
     * G'(t) = lambda*(G(t) - integral_0^t G(t-s)*g(s)ds + S(t))
     * where S(t) = 1 - CDF_g(t) is the survival function of the lifetime distribution.
     * The implicit trapezoidal step resolves to a scalar division at each grid point,
     * giving second-order accuracy and unconditional stability for the decaying modes.</p>
     */
    private void solveVIDE() {
        int N = gridSizeInput.get();
        double h = maxTime / N;

        double[] gridT    = new double[N + 1];
        double[] gridG    = new double[N + 1]; // G(t) = F(t) - 1
        double[] gridGPrime = new double[N + 1]; // G'(t) = F'(t)
        double[] gValues  = new double[N + 1]; // lifetime PDF at grid points
        double[] sValues  = new double[N + 1]; // survival function S(t) = 1 - CDF_g(t)

        org.apache.commons.math3.distribution.RealDistribution dist =
                (org.apache.commons.math3.distribution.RealDistribution) lifetimeDistributionInput.get().getDistribution();

        for (int j = 0; j <= N; j++) {
            gridT[j]   = j * h;
            gValues[j] = dist.density(gridT[j]);
            sValues[j] = 1.0 - dist.cumulativeProbability(gridT[j]);
        }

        // G(0) = 0; G'(0) = lambda*(G(0) - 0 + S(0)) = lambda
        gridG[0]     = 0.0;
        gridGPrime[0] = birthRate;

        // Implicit trapezoidal denominator (constant across steps):
        // derived from the trapezoidal convolution endpoint term 0.5*h*g(0)*G[k+1]
        double g0    = gValues[0];
        double denom = 1.0 - 0.5 * h * birthRate * (1.0 - 0.5 * h * g0);

        for (int k = 0; k < N; k++) {
            // Known interior convolution terms for t_{k+1}:
            // K = h * sum_{j=1}^{k} G[k+1-j] * g[j]
            // (j=0 endpoint contributes 0.5*h*g0*G[k+1] which is the unknown; j=k+1 gives G[0]=0)
            double K = 0.0;
            for (int j = 1; j <= k; j++) {
                K += gridG[k + 1 - j] * gValues[j];
            }
            K *= h;

            // Implicit trapezoidal step: solve for G[k+1]
            // G[k+1] * denom = G[k] + h/2 * G'[k] + h/2 * lambda * (S[k+1] - K)
            gridG[k + 1] = (gridG[k] + 0.5 * h * gridGPrime[k]
                    + 0.5 * h * birthRate * (sValues[k + 1] - K)) / denom;

            // G'[k+1] = lambda*(G[k+1] - Conv_G[k+1] + S[k+1])
            // Conv_G[k+1] = 0.5*h*g0*G[k+1] + K
            double convKp1 = 0.5 * h * g0 * gridG[k + 1] + K;
            gridGPrime[k + 1] = birthRate * (gridG[k + 1] - convKp1 + sValues[k + 1]);
        }

        SplineInterpolator interp = new SplineInterpolator();
        gSpline      = interp.interpolate(gridT, gridG);
        fPrimeSpline = interp.interpolate(gridT, gridGPrime);
    }

    // -------------------------------------------------------------------------
    // Density and CDF
    // -------------------------------------------------------------------------

    @Override
    public double calculateLogNodeAgeDensity(double time) {
        if (lifetimesAreErlang) {
            double[] s       = computeScaledSums(time);
            double maxExp    = s[0], scaledF = s[1], scaledFP = s[2];
            double constTerm = (1.0 - rho) + rho * gammaConst;
            // F_p(t)  = exp(maxExp) * innerFp,   innerFp = constTerm*exp(-maxExp) + rho*scaledF
            // F_p'(t) = rho * exp(maxExp) * scaledFP
            // log(f_p) = log(rho) + log(scaledFP) - maxExp - 2*log(innerFp)
            double innerFp = constTerm * Math.exp(-maxExp) + rho * scaledF;
            return Math.log(rho) + Math.log(scaledFP) - maxExp - 2.0 * Math.log(innerFp);
        } else if (useNumericalSolver) {
            double g       = gSpline.value(time);    // G(t) = F(t) - 1
            double fp      = 1.0 + rho * g;          // fp = (1-rho) + rho*(G+1) = 1 + rho*G
            double fpPrime = rho * fPrimeSpline.value(time);
            return Math.log(fpPrime) - 2.0 * Math.log(fp);
        }
        return 0.0;
    }

    @Override
    public double calculateLogNodeAgeCDF(double time) {
        if (lifetimesAreErlang) {
            double[] s       = computeScaledSums(time);
            double maxExp    = s[0], scaledF = s[1];
            double constTerm = (1.0 - rho) + rho * gammaConst;
            // fp = exp(maxExp) * innerFp
            // log(1 - 1/fp) = log1p(-exp(-maxExp) / innerFp)
            double innerFp = constTerm * Math.exp(-maxExp) + rho * scaledF;
            return Math.log1p(-Math.exp(-maxExp) / innerFp);
        } else if (useNumericalSolver) {
            double g          = gSpline.value(time); // G(t) = F(t) - 1
            double fpMinusOne = rho * g;              // fp - 1 = rho*G, no cancellation
            double fp         = 1.0 + rho * g;
            return Math.log(fpMinusOne) - Math.log(fp);
        }
        return 0.0;
    }

    // -------------------------------------------------------------------------
    // Utilities
    // -------------------------------------------------------------------------

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
