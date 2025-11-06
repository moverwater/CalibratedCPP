package calibratedcpp.model;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
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

    public Double birthRate;
    public Double deathRate;
    public Double diversificationRate;
    public Double reproductiveNumber;
    public Double turnover;
    public Double rho;

    public double logBirthRate;
    public double logDeathRate;
    public double logDiversificationRate;
    public double logRho;

    public double A;
    public double B;

    public boolean isCritical;

    @Override
    public void initAndValidate() {

        birthRate = safeGet(birthRateInput);
        deathRate = safeGet(deathRateInput);
        diversificationRate = safeGet(diversificationRateInput);
        reproductiveNumber = safeGet(reproductiveNumberInput);
        turnover = safeGet(turnoverInput);
        rho = safeGet(rhoInput);

        int specified = 0;
        for (Double i : new Double[]{birthRate, deathRate, diversificationRate, reproductiveNumber, turnover}) {
            if (i != null) specified++;
        }

        if (rho == null) {
            throw new IllegalArgumentException("rho parameter must be specified.");
        }

        if (specified != 2) {
            throw new IllegalArgumentException("Exactly TWO of {birthRate, deathRate, diversificationRate, reproductiveNumber, turnover} must be specified.");
        }

        // disallow repNumber + turnover
        if (reproductiveNumber != null && turnover != null) {
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover together.");
        }

        updateParameters();
    }

    @Override
    public double calculateLogDensity(double time) {
        double logDensity;
        double rt = diversificationRate * time;

        if (isCritical) {
            // Critical case
            logDensity = logRho + logBirthRate - 2 * Math.log1p(A * time);
        } else if (diversificationRate < 0) {
            // Sub-critical case: use stable form with exp(r * t)
            logDensity = logRho + logBirthRate + 2 * logDiversificationRate + rt
                    - 2 * Math.log(Math.abs(A * Math.exp(rt) + B));
        } else {
            // Supercritical case: formula with exp(-r * t)
            logDensity = logRho + logBirthRate + 2 * logDiversificationRate - rt
                    - 2 * Math.log(A + B * Math.exp(-rt));
        }
        return logDensity;
    }

    @Override
    public double calculateLogCDF(double time) {
        double logCDF;

        if (isCritical) {
            // Critical case
            logCDF = logRho + logBirthRate + Math.log(time) - Math.log1p(A * time);
        } else if (diversificationRate < 0) {
            // Sub-critical case
            double exp_rt = Math.exp(diversificationRate * time); // decays, stable

            logCDF = logRho + logBirthRate
                    + Math.log1p(-exp_rt)     // log(1 - exp(r * t)) stable for r<0
                    - Math.log(-A * exp_rt - B);
        } else {
            // Supercritical case
            double exp_neg_rt = Math.exp(-diversificationRate * time);

            logCDF = logRho + logBirthRate
                    + Math.log1p(-exp_neg_rt)  // log(1 - exp(-r * t)) stable for r>=0
                    - Math.log(A + B * exp_neg_rt);
        }
        return logCDF;
    }

    public void updateParameters() {
        birthRate = safeGet(birthRateInput);
        deathRate = safeGet(deathRateInput);
        diversificationRate = safeGet(diversificationRateInput);
        reproductiveNumber = safeGet(reproductiveNumberInput);
        turnover = safeGet(turnoverInput);
        rho = safeGet(rhoInput);

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
            deathRate = diversificationRate / (reproductiveNumber - 1);
            birthRate = deathRate * reproductiveNumber;
        } else if (diversificationRate != null && turnover != null) {
            birthRate = diversificationRate / (1 - turnover);
            deathRate = birthRate * turnover;
        } else {
            throw new IllegalArgumentException("Unsupported parameter combination.");
        }

        A = rho * birthRate;
        B = birthRate * (1 - rho) - deathRate;

        diversificationRate = birthRate - deathRate;

        isCritical = Math.abs(diversificationRate) < 1e-10;

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

    private Double safeGet(Input<RealParameter> input) {
        RealParameter param = input.get();
        return (param != null) ? param.getValue() : null;
    }
}
