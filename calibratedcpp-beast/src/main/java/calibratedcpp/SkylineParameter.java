package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

/**
 * @author Marcus Overwater
 */

@Description("Input for rate parameters in BirthDeathSkylineModel.")
public class SkylineParameter extends CalculationNode {
    public Input<RealVector<? extends Real>> valuesInput =
            new Input<>("values", "Value of the rates specified from root to present.", Input.Validate.REQUIRED);
    public Input<RealVector<NonNegativeReal>> changeTimesInput =
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
            RealVector<?> vals = valuesInput.get();
            for (int i = 0; i < vals.size(); i++) {
                if (vals.get(i) < 0.0) throw new IllegalArgumentException("The rate " + getID() + " must be a non-negative number.");
            }
        }
        if (changeTimesInput.get() != null) {
            if (changeTimesInput.get().size() != valuesInput.get().size() - 1) {
                throw new IllegalArgumentException("Change times of " + this.getID() + " should have dimension equal to the number of rates minus one.");
            }
            RealVector<?> times = changeTimesInput.get();
            for (int i = 0; i < times.size(); i++) {
                double time = times.get(i);
                if (time < 0.0) throw new IllegalArgumentException("The time " + getID() + " must be a non-negative number.");
                if (isRelative && time > 1.0) throw new IllegalArgumentException("When times are relative changeTimes (" + getID() + ") should be less than 1.");
            }
        }
    }
}
