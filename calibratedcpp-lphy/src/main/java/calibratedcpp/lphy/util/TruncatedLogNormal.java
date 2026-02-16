package calibratedcpp.lphy.util;

import lphy.base.distribution.ParametricDistribution;
import lphy.core.model.GenerativeDistribution1D;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.ValueUtils;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Map;
import java.util.TreeMap;

/**
 * A Log-Normal distribution truncated at an upper bound (Right-Truncated).
 * Domain: (0, max]
 * @author Marcus Overwater
 */
public class TruncatedLogNormal extends ParametricDistribution<Double> implements GenerativeDistribution1D<Double> {

    public static final String meanLogParamName = "meanlog";
    public static final String sdLogParamName = "sdlog";
    public static final String maxParamName = "max";

    private Value<Number> M;
    private Value<Number> S;
    private Value<Number> max;

    LogNormalDistribution logNormalDistribution;

    public TruncatedLogNormal(@ParameterInfo(name = meanLogParamName, narrativeName = "mean in log space", description = "the mean of the underlying distribution on the log scale.") Value<Number> M,
                              @ParameterInfo(name = sdLogParamName, narrativeName = "standard deviation in log space", description = "the standard deviation of the underlying distribution on the log scale.") Value<Number> S,
                              @ParameterInfo(name = maxParamName, optional = true, narrativeName = "maximum", description = "the upper bound of the truncation. default is infinity.") Value<Number> max) {
        super();
        this.M = M;
        this.S = S;
        this.max = max;

        constructDistribution(random);
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {
        logNormalDistribution = new LogNormalDistribution(random, ValueUtils.doubleValue(M), ValueUtils.doubleValue(S),
                LogNormalDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
    }

    @GeneratorInfo(name = "TruncatedLogNormal", verbClause = "has", narrativeName = "truncated log-normal distribution",
            category = GeneratorCategory.PRIOR, description = "The log-normal probability distribution truncated to the interval (0, max].")
    public RandomVariable<Double> sample() {
        double maxVal = getMaxVal();

        double pMax = (maxVal == Double.POSITIVE_INFINITY) ? 1.0 : logNormalDistribution.cumulativeProbability(maxVal);

        // 1. Draw u ~ U[0,1]
        double u = random.nextDouble();
        // 2. Scale u to the valid cumulative probability range [0, pMax]
        double p = u * pMax;
        // 3. Map back to x using Inverse CDF
        double result = logNormalDistribution.inverseCumulativeProbability(p);

        return new RandomVariable<>(null, result, this);
    }

    public double logDensity(Double x) {
        double maxVal = getMaxVal();

        if (x <= 0 || x > maxVal) {
            return Double.NEGATIVE_INFINITY;
        }

        double logPdfUnnormalized = logNormalDistribution.logDensity(x);

        double pMax = (maxVal == Double.POSITIVE_INFINITY) ? 1.0 : logNormalDistribution.cumulativeProbability(maxVal);

        return logPdfUnnormalized - Math.log(pMax);
    }

    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(meanLogParamName, M);
            put(sdLogParamName, S);
            if (max != null) put(maxParamName, max);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        switch (paramName) {
            case meanLogParamName:
                M = value;
                break;
            case sdLogParamName:
                S = value;
                break;
            case maxParamName:
                max = value;
                break;
            default:
                throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
        super.setParam(paramName, value);
    }

    private double getMaxVal() {
        return (max != null) ? ValueUtils.doubleValue(max) : Double.POSITIVE_INFINITY;
    }

    public Value<Number> getMeanLog() {
        return M;
    }

    public Value<Number> getSDLog() {
        return S;
    }

    public Double[] getDomainBounds() {
        return new Double[]{0.0, getMaxVal()};
    }
}