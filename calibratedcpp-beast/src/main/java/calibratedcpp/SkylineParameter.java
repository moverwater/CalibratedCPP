package calibratedcpp;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Input for rate parameters in BirthDeathSkylineModel.")
public class SkylineParameter extends BEASTObject {
    public Input<RealParameter> ratesInput =
            new Input<>("rates", "Value of the rates specified from root to present.", Input.Validate.REQUIRED);
    public Input<RealParameter> changeTimesInput =
            new Input<>("changeTimes", "Value of the change times.", Input.Validate.OPTIONAL);
    public Input<Boolean> relativeInput =
            new Input<>("relative", "Boolean whether change times are relative to root height/origin. Default: false", false);
    public Input<Boolean> reverseTimeInput =
            new Input<>("reverseTime", "Boolean: true if change times are specified from present to root. Default: false", false);

    public boolean isRelative;
    public boolean isReverse;

    @Override
    public void initAndValidate() {
        isReverse = reverseTimeInput.get();
        isRelative = relativeInput.get();

        if (ratesInput.get() != null) {
            for (Double rate : ratesInput.get().getValues()){
                if (rate < 0.0) throw new IllegalArgumentException("The rate " + ratesInput.get().getID() + " must be a non-negative number.");
            }
        }
        if (changeTimesInput.get() != null) {
            if (changeTimesInput.get().getDimension() != ratesInput.get().getDimension() - 1) {
                throw new IllegalArgumentException("Change times of " + this.getID() + " should have dimension equal to the number of rates minus one.");
            }
            for (Double time : changeTimesInput.get().getValues()){
                if (time < 0.0) throw new IllegalArgumentException("The time " + changeTimesInput.get().getID() + " must be a non-negative number.");
                if (isRelative && time > 1.0) throw new IllegalArgumentException("When times are relative changeTimes (" + changeTimesInput.get().getID() + ") should be less than 1.");
            }
        }
    }
}
