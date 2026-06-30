package calibratedcpp.lphy.prior;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

public class CalibrationFunction extends DeterministicFunction<Calibration> {

    public static final String taxaParamName = "taxa";
    public static final String upperParamName = "upper";
    public static final String lowerParamName = "lower";

    public CalibrationFunction(
            @ParameterInfo(name = taxaParamName, description = "the taxa defining the MRCA node to calibrate") Value<String[]> taxa,
            @ParameterInfo(name = upperParamName, description = "the upper bound of the calibration age") Value<Number> upper,
            @ParameterInfo(name = lowerParamName, description = "the lower bound of the calibration age") Value<Number> lower) {
        setParam(taxaParamName, taxa);
        setParam(upperParamName, upper);
        setParam(lowerParamName, lower);
    }

    @GeneratorInfo(name = "calibration",
            description = "Creates a calibration constraint for the MRCA of the given taxa, with the specified age bounds.")
    @Override
    public Value<Calibration> apply() {
        String[] taxa = ((Value<String[]>) getParams().get(taxaParamName)).value();
        double upper = ((Value<Number>) getParams().get(upperParamName)).value().doubleValue();
        double lower = ((Value<Number>) getParams().get(lowerParamName)).value().doubleValue();
        return new Value<>(null, new Calibration(taxa, upper, lower), this);
    }
}
