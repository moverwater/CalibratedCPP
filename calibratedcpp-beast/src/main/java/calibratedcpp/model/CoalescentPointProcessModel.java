package calibratedcpp.model;

import beast.base.core.Description;
import beast.base.inference.CalculationNode;

/**
 * @author Marcus Overwater
 */

@Description("Abstract class for the distribution of node ages in a CPP")
public abstract class CoalescentPointProcessModel extends CalculationNode {
    public abstract double calculateLogDensity(double time);

    public abstract double calculateLogCDF(double time);
}