package calibratedcpp.model;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * Implements the node age distribution under a constant-rate birth-death process
 * for the Coalescent Point Process (CPP) representation.
 *
 * <p>This model defines how node ages are distributed in a phylogenetic tree
 * generated under a constant birth-death process, parameterized by combinations of
 * rates such as birth (λ), death (μ), diversification (λ - μ), reproductive number (λ / μ),
 * turnover (μ / λ), and sampling probability (ρ).</p>
 *
 * <p>Exactly two parameters among {birthRate, deathRate, diversificationRate,
 * reproductiveNumber, turnover} must be provided, along with rho (sampling probability).
 * The remaining parameters are derived automatically.</p>
 *
 * <p>Based on analytical derivations of CPP distributions for birth-death models.</p>
 *
 * @author Marcus Overwater
 */
@Description("Node age distribution for the CPP representation of the birth-death process")
public class BirthDeathModel extends CoalescentPointProcessModel {

    public Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "The birth rate (lambda)", (RealParameter) null);

    public Input<RealParameter> deathRateInput =
            new Input<>("deathRate", "The death rate (mu)", (RealParameter) null);

    public Input<RealParameter> diversificationRateInput =
            new Input<>("diversificationRate", "Diversification rate (lambda - mu)", (RealParameter) null);

    public Input<RealParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "Reproductive number (lambda / mu)", (RealParameter) null);

    public Input<RealParameter> turnoverInput =
            new Input<>("turnover", "Turnover (mu / lambda)", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho", "Probability with which each individual in the population is sampled.", (RealParameter) null);

    // Numeric versions of model parameters after extraction from inputs
    public Double birthRate;
    public Double deathRate;
    public Double diversificationRate;
    public Double reproductiveNumber;
    public Double turnover;
    public Double rho;

    // Cached logarithmic values for efficiency
    public double logBirthRate;
    public double logDeathRate;
    public double logDiversificationRate;
    public double logRho;

    // Derived constants for density and CDF calculations
    public double A;
    public double B;

    // Indicates whether diversification rate ≈ 0 (critical case)
    public boolean isCritical;

    @Override
    public void initAndValidate() {

        // Extract initial parameter values safely
        birthRate = safeGet(birthRateInput);
        deathRate = safeGet(deathRateInput);
        diversificationRate = safeGet(diversificationRateInput);
        reproductiveNumber = safeGet(reproductiveNumberInput);
        turnover = safeGet(turnoverInput);
        rho = safeGet(rhoInput);

        // Count how many of the five possible parameters were provided
        int specified = 0;
        StringBuilder whichSpecified = new StringBuilder();

        Input<?>[] rateInputs = {birthRateInput, deathRateInput, reproductiveNumberInput, diversificationRateInput, turnoverInput};

        for (Input<?> input : rateInputs) {
            if (input.get() != null) {
                specified++;
                if (!whichSpecified.isEmpty()) whichSpecified.append(", ");
                whichSpecified.append(input.getName());
            }
        }

        // rho is mandatory
        if (rho == null) {
            throw new IllegalArgumentException("rho parameter must be specified.");
        }

        if (specified != 2) {
            throw new IllegalArgumentException("Exactly TWO of {birthRate, deathRate, reproductiveNumber, diversificationRate, turnover} must be specified. " +
                    "These " + specified + " were specified: (" + whichSpecified + ")");
        }

        // Disallow incompatible combinations
        if (reproductiveNumber != null && turnover != null) {
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover together.");
        }

        // Compute all derived parameters
        updateParameters();
    }

    /**
     * Computes the log of the probability density function (PDF)
     * for a given node age {@code time}.
     *
     * @param time node age
     * @return log-density value
     */
    @Override
    public double calculateLogDensity(double time) {
        double logDensity;
        double rt = diversificationRate * time;

        if (isCritical) {
            // Case when diversification rate ≈ 0
            logDensity = logRho + logBirthRate - 2 * Math.log1p(A * time);
        } else if (diversificationRate < 0) {
            // Subcritical case: r < 0
            logDensity = logRho + logBirthRate + 2 * logDiversificationRate + rt
                    - 2 * Math.log(Math.abs(A * Math.exp(rt) + B));
        } else {
            // Supercritical case: r > 0
            logDensity = logRho + logBirthRate + 2 * logDiversificationRate - rt
                    - 2 * Math.log(A + B * Math.exp(-rt));
        }
        return logDensity;
    }

    /**
     * Computes the log of the cumulative distribution function (CDF)
     * for a given node age {@code time}.
     *
     * @param time node age
     * @return log-CDF value
     */
    @Override
    public double calculateLogCDF(double time) {
        double logCDF;

        if (isCritical) {
            // Critical case: r ≈ 0
            logCDF = logRho + logBirthRate + Math.log(time) - Math.log1p(A * time);
        } else if (diversificationRate < 0) {
            // Subcritical case: r < 0
            double exp_rt = Math.exp(diversificationRate * time); // decaying exponential
            logCDF = logRho + logBirthRate
                    + Math.log1p(-exp_rt)
                    - Math.log(-A * exp_rt - B);
        } else {
            // Supercritical case: r > 0
            double exp_neg_rt = Math.exp(-diversificationRate * time);
            logCDF = logRho + logBirthRate
                    + Math.log1p(-exp_neg_rt)
                    - Math.log(A + B * exp_neg_rt);
        }
        return logCDF;
    }

    /**
     * Updates derived parameters based on which two input parameters were specified.
     *
     * <p>This method determines the remaining model parameters using algebraic relationships
     * among λ (birth), μ (death), diversification (λ - μ), reproductive number (λ / μ),
     * and turnover (μ / λ).</p>
     *
     * <p>It also computes derived constants (A, B), criticality flag,
     * and cached logarithmic values for efficiency.</p>
     */
    public void updateParameters() {
        // Retrieve current parameter values
        birthRate = safeGet(birthRateInput);
        deathRate = safeGet(deathRateInput);
        diversificationRate = safeGet(diversificationRateInput);
        reproductiveNumber = safeGet(reproductiveNumberInput);
        turnover = safeGet(turnoverInput);
        rho = safeGet(rhoInput);

        // Derive missing parameters depending on which are defined
        if (birthRate != null && deathRate != null) {
            diversificationRate = birthRate - deathRate;
        } else if (birthRate != null && diversificationRate != null) {
            deathRate = birthRate - diversificationRate;
        } else if (deathRate != null && diversificationRate != null) {
            birthRate = deathRate + diversificationRate;
        } else if (birthRate != null && reproductiveNumber != null) {
            deathRate = birthRate / reproductiveNumber;
        } else if (deathRate != null && reproductiveNumber != null) {
            birthRate = deathRate * reproductiveNumber;
        } else if (birthRate != null && turnover != null) {
            deathRate = birthRate * turnover;
        } else if (deathRate != null && turnover != null) {
            birthRate = deathRate / turnover;
        } else if (diversificationRate != null && reproductiveNumber != null) {
            if (Math.abs(reproductiveNumber - 1.0) < 1e-10) {
                throw new IllegalArgumentException("Reproductive number cannot be exactly 1.0 when deriving rates from diversification.");
            }
            deathRate = diversificationRate / (reproductiveNumber - 1);
            birthRate = deathRate * reproductiveNumber;
        } else if (diversificationRate != null && turnover != null) {
            if (Math.abs(1.0 - turnover) < 1e-10) {
                throw new IllegalArgumentException("Turnover cannot be exactly 1.0 when deriving rates from diversification.");
            }
            birthRate = diversificationRate / (1 - turnover);
            deathRate = birthRate * turnover;
        } else {
            throw new IllegalArgumentException("Unsupported parameter combination.");
        }

        // Precompute constants for PDF/CDF
        A = rho * birthRate;
        B = birthRate * (1 - rho) - deathRate;

        diversificationRate = birthRate - deathRate;

        // Check if diversification ≈ 0
        isCritical = Math.abs(diversificationRate) < 1e-10;

        // Cache logarithmic values for performance
        logBirthRate = Math.log(birthRate);
        logDeathRate = Math.log(deathRate);
        logRho = Math.log(rho);
        logDiversificationRate = Math.log(Math.abs(diversificationRate));
    }

    @Override
    protected boolean requiresRecalculation() {
        updateParameters();
        return true;
    }

    /**
     * Safely retrieves a numeric value from an {@link Input<RealParameter>}.
     *
     * @param input the RealParameter input
     * @return the numeric value, or {@code null} if input is undefined
     */
    private Double safeGet(Input<RealParameter> input) {
        RealParameter param = input.get();
        return (param != null) ? param.getValue() : null;
    }

    @Override
    protected void restore() {
        updateParameters();
        super.restore();
    }
}
