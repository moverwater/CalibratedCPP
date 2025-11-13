package calibratedcpp.model;

import beast.base.core.Description;
import beast.base.inference.CalculationNode;

/**
 * Abstract base class representing the distribution of node ages
 * in a Coalescent Point Process (CPP).
 *
 * <p>This class defines the core interface for any model that provides
 * a probability distribution over node ages in a CPP framework.
 * Specific models must extend this class
 * and implement the abstract methods for computing the log-density
 * and log-cumulative distribution functions.</p>
 *
 * <p>The CPP framework provides a statistical representation of coalescent
 * events in evolutionary trees, where node ages are treated as random variables
 * drawn from an underlying generative process (e.g., birth-death, Yule, etc.).</p>
 *
 * @author Marcus Overwater
 */
@Description("Abstract class for the distribution of node ages in a Coalescent Point Process (CPP)")
public abstract class CoalescentPointProcessModel extends CalculationNode {

    /**
     * Calculates the logarithm of the probability density function (PDF)
     * for the age of a node in the CPP.
     *
     * <p>The returned value represents {@code log(q(t))}, where {@code q(t)} is
     * the density of observing a node at age {@code t}} under the model.</p>
     *
     * @param time the node age (time since origin or present)
     * @return the log of the probability density for the given {@code time}
     */
    public abstract double calculateLogDensity(double time);

    /**
     * Calculates the logarithm of the cumulative distribution function (CDF)
     * for the age of a node in the CPP.
     *
     * <p>The returned value represents {@code log(Q(t))}, where {@code Q(t)} is
     * the probability that a node age is less than or equal to {@code t}}.</p>
     *
     * @param time the node age (time since origin or present)
     * @return the log of the cumulative probability for the given {@code time}
     */
    public abstract double calculateLogCDF(double time);
}
