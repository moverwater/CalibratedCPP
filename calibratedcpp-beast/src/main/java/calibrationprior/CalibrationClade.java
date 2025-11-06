package calibrationprior;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.parameter.RealParameter;

/**
 * @author Marcus Overwater
 */

@Description("A calibration clade is a taxon set in a monophyletic clade and an upper and lower bound on the age of the clade")
public class CalibrationClade extends BEASTObject {
    Input<TaxonSet> taxa =
            new Input<>("taxa","the set of taxa in the clade", Input.Validate.REQUIRED);
    Input<RealParameter> upperAgeInput =
            new Input<>("upperAge","the soft upper bound on the age of the clade", Input.Validate.REQUIRED);
    Input<RealParameter> lowerAgeInput =
            new Input<>("lowerAge","the soft lower bound on the age of the clade", Input.Validate.REQUIRED);
    Input<RealParameter> pCoverageInput =
            new Input<>("confidenceLevel","the amount of probability mass in the bounds" +
            "default value (0.9)", new RealParameter("0.9"));

    @Override
    public void initAndValidate() {
        double p = pCoverageInput.get().getValue();
        if ((p < 0.0) || (p >= 1.0)){
            throw new IllegalArgumentException("confidenceLevel (" + p + ") should be between 0.0 and 1.0");
        }

        double t_lo = lowerAgeInput.get().getValue();
        double t_hi = upperAgeInput.get().getValue();

        if (t_hi < t_lo) {
            throw new IllegalArgumentException("lowerAge (" + t_lo + ") should be less than upperAge (" + t_hi + ")");
        }
        if ((t_lo < 0.0) || (t_hi < 0.0)) {
            throw new IllegalArgumentException("lowerAge should be greater than 0.0: lowerAge=(" + t_lo + "), upperAge=(" + t_hi + ")" );
        }
    }

    public TaxonSet getTaxa() {
        return taxa.get();
    }

    public double getUpperAge() {
        return upperAgeInput.get().getValue();
    }

    public double getLowerAge() {
        return lowerAgeInput.get().getValue();
    }

    public double getPCoverage() {
        return pCoverageInput.get().getValue();
    }
}
