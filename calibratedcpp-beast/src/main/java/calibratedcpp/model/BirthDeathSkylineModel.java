package calibratedcpp.model;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Node age distribution for the CPP representation of the birth-death process with piecewise constant rates")
public class BirthDeathSkylineModel extends CoalescentPointProcessModel {
    public Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "the birth rate", (RealParameter) null);

    public Input<RealParameter> deathRateInput =
            new Input<>("deathRate", "the death rate", (RealParameter) null);

    public Input<RealParameter> reproductiveNumberInput =
            new Input<>("reproductiveNumber", "the reproductive number, birthRate / deathRate", (RealParameter) null);

    public Input<RealParameter> diversificationRateInput =
            new Input<>("diversificationRate", "the diversification rate, birthRate - deathRate", (RealParameter) null);

    public Input<RealParameter> turnoverInput =
            new Input<>("turnover", "deathRate / birthRate", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho", "the probability with which each individual in the total population is sampled", (RealParameter) null);

    public Input<RealParameter> birthRateChangeTimesInput =
            new Input<>("birthRateChangeTimes", "the birth rate change times", (RealParameter) null);

    public Input<RealParameter> deathRateChangeTimesInput =
            new Input<>("deathRateChangeTimes", "the death rate change times", (RealParameter) null);

    protected double birthRate;
    protected double deathRate;
    protected double diversificationRate;
    protected double reproductiveNumber;
    protected double turnover;
    protected double rho;
    protected double birthRateChangeTimes;
    protected double deathRateChangeTimes;

    @Override
    public void initAndValidate() {
        birthRate = birthRateInput.get().getValue();
        deathRate = deathRateInput.get().getValue();
        reproductiveNumber = reproductiveNumberInput.get().getValue();
        diversificationRate = diversificationRateInput.get().getValue();
        turnover = turnoverInput.get().getValue();
        rho = rhoInput.get().getValue();
        birthRateChangeTimes = birthRateChangeTimesInput.get().getValue();
        deathRateChangeTimes = deathRateChangeTimesInput.get().getValue();
    }

    @Override
    public double calculateLogDensity(double time) {
        return 0;
    }

    @Override
    public double calculateLogCDF(double time) {
        return 0;
    }
}