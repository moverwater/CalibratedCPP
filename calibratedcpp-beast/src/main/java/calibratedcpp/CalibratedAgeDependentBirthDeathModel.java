package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Gamma;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegratorFactory;
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

    public Input<ParametricDistribution> lifetimeDistributionInput = new Input<>("lifetimeDistribution",
            "Distribution of the lifetime of an individual.");
    public Input<RealParameter> rhoInput = new Input<>("rho", "Extant sampling probability.");
    public Input<RealParameter> birthRateInput = new Input<>("birthRate", "The birth rate.");
    public Input<Integer> gridSizeInput = new Input<>("gridSize",
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
    protected PolynomialSplineFunction lifetimePdfSpline;  // g(s): lifetime PDF
    protected PolynomialSplineFunction survivalSpline;     // S(t) = 1 - CDF_g(t)

    private static final GaussIntegratorFactory GAUSS_FACTORY = new GaussIntegratorFactory();

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
//        if (Math.abs(shapeParam - n) > 1e-10 || n < 1) {
//            throw new IllegalArgumentException(
//                    "Shape parameter must be a positive integer for Erlang distribution, got " + shapeParam);
//        }
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
     * Solves the VIDE G'(t) = lambda*(G(t) - integral_0^t G(t-s)*g(s)ds + S(t))
     * on [0, maxTime] and builds cubic splines for G, the lifetime PDF, and the survival function.
     *
     * <p>Uses Richardson extrapolation over two implicit-trapezoidal solves (step h and step 2h)
     * to cancel the O(h^2) leading error, giving O(h^4) global accuracy in G(t) at ~25% extra cost.
     * The combined G values are at the coarse-grid points (N/2 + 1 points), which is sufficient
     * resolution for the cubic spline.</p>
     */
    private void solveVIDE() {
        int N = gridSizeInput.get();
        if (N % 2 != 0) N++;  // Richardson requires an even step count

        ParametricDistribution dist = lifetimeDistributionInput.get();

        double[] gridG_fine   = computeGridG(N,     dist);
        double[] gridG_coarse = computeGridG(N / 2, dist);

        // Richardson extrapolation: G_4th = (4*G_h - G_2h) / 3 at coarse-grid points.
        // Splines for the lifetime PDF and survival function are exact (no Richardson needed).
        int M = N / 2;
        double hc = maxTime / M;
        double[] coarseT  = new double[M + 1];
        double[] gridG    = new double[M + 1];
        double[] gValues  = new double[M + 1];
        double[] sValues  = new double[M + 1];

        try {
            for (int j = 0; j <= M; j++) {
                coarseT[j] = j * hc;
                gridG[j]   = (4.0 * gridG_fine[2 * j] - gridG_coarse[j]) / 3.0;
                gValues[j] = dist.density(coarseT[j]);
                sValues[j] = 1.0 - dist.cumulativeProbability(coarseT[j]);
            }
        } catch (org.apache.commons.math.MathException e) {
            throw new RuntimeException("Failed to evaluate lifetime distribution", e);
        }

        SplineInterpolator interp = new SplineInterpolator();
        gSpline           = interp.interpolate(coarseT, gridG);
        lifetimePdfSpline = interp.interpolate(coarseT, gValues);
        survivalSpline    = interp.interpolate(coarseT, sValues);
    }

    /**
     * Runs the implicit trapezoidal VIDE solver on a uniform grid of N steps over [0, maxTime].
     * Returns G(t) = F(t) - 1 at all N+1 grid points.
     */
    private double[] computeGridG(int N, ParametricDistribution dist) {
        double h = maxTime / N;
        double[] gridG      = new double[N + 1];
        double[] gridGPrime = new double[N + 1];
        double[] gValues    = new double[N + 1];
        double[] sValues    = new double[N + 1];

        try {
            for (int j = 0; j <= N; j++) {
                double t   = j * h;
                gValues[j] = dist.density(t);
                sValues[j] = 1.0 - dist.cumulativeProbability(t);
            }
        } catch (org.apache.commons.math.MathException e) {
            throw new RuntimeException("Failed to evaluate lifetime distribution at a grid point", e);
        }

        gridG[0]      = 0.0;
        gridGPrime[0] = birthRate;
        double g0    = gValues[0];
        double denom = 1.0 - 0.5 * h * birthRate * (1.0 - 0.5 * h * g0);

        for (int k = 0; k < N; k++) {
            double K = 0.0;
            for (int j = 1; j <= k; j++) {
                K += gridG[k + 1 - j] * gValues[j];
            }
            K *= h;

            gridG[k + 1] = (gridG[k] + 0.5 * h * gridGPrime[k]
                    + 0.5 * h * birthRate * (sValues[k + 1] - K)) / denom;

            double convKp1    = 0.5 * h * g0 * gridG[k + 1] + K;
            gridGPrime[k + 1] = birthRate * (gridG[k + 1] - convKp1 + sValues[k + 1]);
        }
        return gridG;
    }

    /**
     * Evaluates G'(t) = lambda*(G(t) - integral_0^t G(t-s)*g(s)ds + S(t)) on demand.
     *
     * <p>Rather than interpolating a stored G'(t) spline (which accumulates O(h^2) absolute
     * error and gives large relative error when G'(t) is exponentially small in the subcritical
     * case), we re-evaluate the VIDE formula directly using 32-point Gauss-Legendre quadrature
     * on the G spline. The O(h^2) equilibrium error in G(t) cancels with the corresponding
     * error in the convolution integral, leaving a remainder of order O(h^2 * S(t)) which
     * decays exponentially — giving stable relative accuracy at all t.</p>
     */
    private double evaluateGPrime(double t) {
        if (t <= 0.0) return birthRate;
        GaussIntegrator gl = GAUSS_FACTORY.legendre(32, 0.0, t);
        double integral = gl.integrate(s -> gSpline.value(t - s) * lifetimePdfSpline.value(s));
        return birthRate * (gSpline.value(t) - integral + survivalSpline.value(t));
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
            double g       = gSpline.value(time);
            double fp      = 1.0 + rho * g;
            double fpPrime = rho * evaluateGPrime(time);
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
