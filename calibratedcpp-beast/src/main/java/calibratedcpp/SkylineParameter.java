package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("Input for rate parameters in BirthDeathSkylineModel.")
public class SkylineParameter extends CalculationNode {
    public Input<RealParameter> valuesInput =
            new Input<>("values", "Value of the rates specified from root to present.", Input.Validate.REQUIRED);
    public Input<RealParameter> changeTimesInput =
            new Input<>("changeTimes", "Value of the change times.", Input.Validate.OPTIONAL);
    public Input<Boolean> timesAreRelativeInput =
            new Input<>("timesAreRelative", "Boolean whether change times are relative to root height/origin. Default: false", false);
    public Input<Boolean> timesAreAgesInput =
            new Input<>("timesAreAges", "Boolean: true if change times are specified from present to root. Default: false", false);

    public boolean isRelative;
    public boolean isReverse;

    @Override
    public void initAndValidate() {
        isReverse = timesAreAgesInput.get();
        isRelative = timesAreRelativeInput.get();

        if (valuesInput.get() != null) {
            for (Double rate : valuesInput.get().getValues()){
                if (rate < 0.0) throw new IllegalArgumentException("The rate " + valuesInput.get().getID() + " must be a non-negative number.");
            }
        }
        if (changeTimesInput.get() != null) {
            if (changeTimesInput.get().getDimension() != valuesInput.get().getDimension() - 1) {
                throw new IllegalArgumentException("Change times of " + this.getID() + " should have dimension equal to the number of rates minus one.");
            }
            for (Double time : changeTimesInput.get().getValues()){
                if (time < 0.0) throw new IllegalArgumentException("The time " + changeTimesInput.get().getID() + " must be a non-negative number.");
                if (isRelative && time > 1.0) throw new IllegalArgumentException("When times are relative changeTimes (" + changeTimesInput.get().getID() + ") should be less than 1.");
            }
        }
    }
}
