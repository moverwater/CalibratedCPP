package calibratedcpp.lphy.prior;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

public class toCalibrationArray extends DeterministicFunction<CalibrationArray> {
    public final String calibration1Name = "calibration1";
    public final String calibration2Name = "calibration2";

    public toCalibrationArray(@ParameterInfo(name = calibration1Name,description = "calibrations array")Value<Calibration> calibration1,
                              @ParameterInfo(name = calibration2Name, description = "") Value<Calibration> calibration2) {
        setParam(calibration1Name, calibration1);
        setParam(calibration2Name, calibration2);
    }

    @GeneratorInfo(name = "toArray", description = "")
    @Override
    public Value<CalibrationArray> apply() {
        Value<Calibration> c1 = getParams().get(calibration1Name);
        Value<Calibration> c2 = getParams().get(calibration2Name);
        Calibration calibration1 = c1.value();
        Calibration calibration2 = c2.value();

        Calibration[] calibrationsValue = new Calibration[] {calibration1, calibration2};
        CalibrationArray calibrationArray = new CalibrationArray(calibrationsValue);
        return new Value<>("", calibrationArray, this);
    }
}
