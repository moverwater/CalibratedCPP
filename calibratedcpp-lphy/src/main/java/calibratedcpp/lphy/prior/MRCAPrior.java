package calibratedcpp.lphy.prior;

import lphy.base.distribution.ParametricDistribution;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;
import java.util.TreeMap;

public class MRCAPrior implements GenerativeDistribution<Calibration> {

    public static final String distributionParamName = "distribution";
    public static final String taxaParamName = "taxa";

    private Value<ParametricDistribution<Double>> distr;
    private Value<String[]> taxonSet;

    public MRCAPrior(@ParameterInfo(name = distributionParamName, narrativeName = "distribution of TMRCA", description = "The distribution of the tmrca of the clade.") Value<ParametricDistribution<Double>> distr,
                     @ParameterInfo(name = taxaParamName, narrativeName = "set of taxa", description = "The set of taxa.") Value<String[]> taxonSet) {
        this.distr = distr;
        this.taxonSet = taxonSet;
    }

    @Override
    public RandomVariable<Calibration> sample() {
        // 1. Get the underlying distribution object
        ParametricDistribution<Double> distribution = distr.value();

        // 2. Sample a random age (Double) from that distribution
        Double sampledAge = distribution.sample().value();

        // 3. Create the Calibration object using the provided taxa and sampled age
        Calibration calibration = new Calibration(taxonSet.value(), sampledAge);

        // 4. Return as a RandomVariable<Calibration>
        return new RandomVariable<>(null, calibration, this);
    }

    @Override
    public double logDensity(Calibration value) {
        if (value.getAge() == null) {
            return Double.NEGATIVE_INFINITY;
        }
        return distr.value().logDensity(value.getAge());
    }

    @Override
    public Map<String, Value> getParams() {
        // Use TreeMap to ensure consistent parameter ordering
        return new TreeMap<>() {{
            put(distributionParamName, distr);
            put(taxaParamName, taxonSet);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(distributionParamName)) {
            this.distr = value;
        } else if (paramName.equals(taxaParamName)) {
            this.taxonSet = value;
        } else {
            throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
    }
}