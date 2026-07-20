package calibratedcpp.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;

/**
 * @author Marcus Overwater
 */
@Description("Proposal operator for skyline change times which preserves their ordering.")
public class ChangeTimeOperator extends Operator {

    public Input<RealVectorParam<? extends Real>> changeTimesInput = new Input<>("changeTimes",
            "Change times to operate on; must be strictly increasing.", Input.Validate.REQUIRED);

    public Input<Double> windowSizeInput = new Input<>("windowSize",
            "Size of the window from which to draw a shift for the first change time.", 1.0);

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Used to scale the gaps between adjacent times. The scale factor is drawn from " +
                    "between scaleFactor and 1/scaleFactor.", 0.75);

    private RealVectorParam<? extends Real> changeTimes;
    private double lower, upper;

    @Override
    public void initAndValidate() {
        changeTimes = changeTimesInput.get();

        // Bounds come from the parameter's domain rather than from operator inputs. In beast-base
        // 2.8 the spec parameter classes no longer carry their own lower/upper (those inputs are
        // disabled), so the domain is where a hard support lives; NonNegativeReal gives [0, +inf).
        // These checks are only a short-circuit that avoids a wasted likelihood evaluation: a
        // proposal outside the prior's support is rejected by the MH ratio regardless, which is why
        // an unbounded domain here is safe.
        Real domain = changeTimes.getDomain();
        lower = domain.getLower();
        upper = domain.getUpper();

        for (int i = 1; i < changeTimes.size(); i++) {
            if (changeTimes.get(i) <= changeTimes.get(i - 1))
                throw new IllegalArgumentException("ChangeTimeOperator can only be applied to a " +
                        "strictly increasing sequence, but " + changeTimes.getID() + " entry " + i +
                        " (" + changeTimes.get(i) + ") is not greater than entry " + (i - 1) +
                        " (" + changeTimes.get(i - 1) + ").");
        }
    }

    @Override
    public double proposal() {
        int idx = Randomizer.nextInt(changeTimes.size());
        return idx == 0 ? shiftProposal() : gapProposal(idx);
    }

    /** Shifts the whole sequence by a uniform window draw, preserving all gaps. */
    private double shiftProposal() {
        double delta = (Randomizer.nextDouble() - 0.5) * windowSizeInput.get();

        if (changeTimes.get(0) + delta < lower
                || changeTimes.get(changeTimes.size() - 1) + delta > upper)
            return Double.NEGATIVE_INFINITY;

        for (int i = 0; i < changeTimes.size(); i++)
            changeTimes.set(i, changeTimes.get(i) + delta);

        return 0.0;
    }

    /** Scales the gap between idx-1 and idx, shifting idx and everything after it. */
    private double gapProposal(int idx) {
        double minf = Math.min(scaleFactorInput.get(), 1.0 / scaleFactorInput.get());
        double f = minf + Randomizer.nextDouble() * (1 / minf - minf);

        double shift = (f - 1) * (changeTimes.get(idx) - changeTimes.get(idx - 1));

        // Checked after the shift, not before: BDMM-Prime tests the pre-shift value here, which
        // lets an out-of-bounds proposal through (harmlessly, since the prior then rejects it).
        if (changeTimes.get(changeTimes.size() - 1) + shift > upper)
            return Double.NEGATIVE_INFINITY;

        for (int i = idx; i < changeTimes.size(); i++)
            changeTimes.set(i, changeTimes.get(i) + shift);

        return -Math.log(f);
    }
}
