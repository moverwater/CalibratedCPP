package calibratedcpp.lphy.prior;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import static calibratedcpp.lphy.prior.ConditionedPriorUtils.extractCalibrations;

/**
 * Collects any number of individual {@link Calibration} values (e.g. from {@link UniformMRCA})
 * into one {@link CalibrationArray}, using the same array-literal syntax as
 * {@code ConditionedMRCAPrior(calibrations=[cal1, cal2, ...])}:
 * {@code toArray(calibrations=[cal1, cal2, ..., calN])}. Scales to any number of calibrations
 * without pairwise nesting.
 */
public class toCalibrationArray extends DeterministicFunction<CalibrationArray> {
    public static final String calibrationsParamName = "calibrations";

    public toCalibrationArray(@ParameterInfo(name = calibrationsParamName, description = "array of individual calibration constraints") Value<?> calibrations) {
        if (calibrations == null)
            throw new IllegalArgumentException("calibrations must be provided");
        setParam(calibrationsParamName, calibrations);
    }

    @GeneratorInfo(name = "toArray", description = "Collects individual calibrations into one CalibrationArray.")
    @Override
    public Value<CalibrationArray> apply() {
        Value<?> calibrations = getParams().get(calibrationsParamName);
        CalibrationArray calibrationArray = new CalibrationArray(extractCalibrations(calibrations));
        return new Value<>("", calibrationArray, this);
    }
}
