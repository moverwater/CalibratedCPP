package calibratedcpp.lphy.prior;

import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;
import java.util.TreeMap;

/*
    This generator is just for validation of tree simulator, not usable when doing actual experiments
 */
public class MRCAPrior implements GenerativeDistribution<Calibration> {

    public static final String ageName = "age";
    public static final String taxaParamName = "taxa";

    private Value<Double> age;
    private Value<String[]> taxonSet;

    public MRCAPrior(@ParameterInfo(name = ageName, narrativeName = "distribution of TMRCA", description = "The distribution of the tmrca of the clade.") Value<Double> age,
                     @ParameterInfo(name = taxaParamName, narrativeName = "set of taxa", description = "The set of taxa.") Value<String[]> taxonSet) {
        this.age = age;
        this.taxonSet = taxonSet;
    }

    @GeneratorInfo(name = "MRCAPrior", description = "")
    @Override
    public RandomVariable<Calibration> sample() {
        double age = getAge().value();

        // 3. Create the Calibration object using the provided taxa and sampled age
        Calibration calibration = new Calibration(taxonSet.value(), age);

        // 4. Return as a RandomVariable<Calibration>
        return new RandomVariable<>(null, calibration, this);
    }

    @Override
    public Map<String, Value> getParams() {
        // Use TreeMap to ensure consistent parameter ordering
        return new TreeMap<>() {{
            put(ageName, age);
            put(taxaParamName, taxonSet);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(ageName)) {
            this.age = value;
        } else if (paramName.equals(taxaParamName)) {
            this.taxonSet = value;
        } else {
            throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
    }

    public Value<Double> getAge(){
        return getParams().get(ageName);
    }
}