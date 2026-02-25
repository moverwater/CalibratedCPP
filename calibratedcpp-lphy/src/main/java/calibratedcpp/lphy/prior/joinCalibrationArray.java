package calibratedcpp.lphy.prior;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

/*
    This class is just for type conversion from Calibration to CalibrationArray, non-usable in real experiments
 */
public class joinCalibrationArray extends DeterministicFunction<CalibrationArray> {
    public final String calibration1Name = "calibration1";
    public final String calibration2Name = "calibration2";
    public final String calibration3Name = "singleCalibration";

    public joinCalibrationArray(@ParameterInfo(name = calibration1Name,description = "calibrations array")Value<CalibrationArray> calibration1,
                                @ParameterInfo(name = calibration2Name, description = "", optional = true) Value<CalibrationArray> calibration2,
                                @ParameterInfo(name = calibration3Name, description = "", optional = true) Value<Calibration> singleCalibration) {
        setParam(calibration1Name, calibration1);
        if (calibration2 != null) {
            setParam(calibration2Name, calibration2);
        }
        if(singleCalibration != null) {
            setParam(calibration3Name, singleCalibration);
        }
        if (calibration2 == null && singleCalibration == null) {
            throw new IllegalArgumentException("Both calibrations and singleCalibration are null, should give at least 1");
        }
    }

    @GeneratorInfo(name = "joinArray", description = "")
    @Override
    public Value<CalibrationArray> apply() {
        Value<CalibrationArray> c1 = getParams().get(calibration1Name);
        Value<CalibrationArray> c2 = getParams().get(calibration2Name);
        Value<Calibration> c3 = getParams().get(calibration3Name);
        Calibration[] calibration1 = c1.value().array;
        CalibrationArray calibrationArray;
        if (c2!= null) {
            Calibration[] calibration2 = c2.value().array;
            Calibration[] calibrationsValue = new Calibration[calibration1.length + calibration2.length];
            for (int i = 0; i < calibration1.length; i++) {
                calibrationsValue[i] = calibration1[i];
            }
            for (int i = 0; i < calibration2.length; i++) {
                calibrationsValue[i + calibration1.length] = calibration2[i];
            }
            calibrationArray = new CalibrationArray(calibrationsValue);
        } else {
            calibrationArray = new CalibrationArray(calibration1);
        }
        if (c3 != null) {
            calibrationArray.addCalibration(c3.value());
        }

        return new Value<>("", calibrationArray, this);
    }
}
