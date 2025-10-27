package calibratedcpp.model;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("This gives the node age distribution and density of an age dependent binary branching process" +
        "where individuals share some time dependent lifetime distribution and give birth at a piecewise constant rate.")
public class AgeDependentModel extends CoalescentPointProcessModel {
    public Input<Distribution> lifetimeDistributionInput =
            new Input<>("lifetimeDistribution", "distribution of the lifetime of an individual", (Distribution) null);

    public Input<RealParameter> birthRateInput =
            new Input<>("birthRate", "the rate at which individuals give birth", (RealParameter) null);

    protected double birthRate;
    protected Distribution lifetimeDistribution;

    @Override
    public void initAndValidate() {
        lifetimeDistribution = lifetimeDistributionInput.get();
        birthRate = birthRateInput.get().getValue();
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
