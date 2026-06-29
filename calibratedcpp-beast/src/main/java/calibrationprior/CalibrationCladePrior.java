package calibrationprior;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;

/**
 * @author Marcus Overwater
 */

@Description("A calibration clade with a soft upper and lower bound on the age of the clade MRCA.")
public class CalibrationCladePrior extends BEASTObject {
    public Input<TaxonSet> taxaInput =
            new Input<>("taxa", "The set of taxa whose MRCA age is calibrated.", Input.Validate.REQUIRED);
    public Input<RealScalar<NonNegativeReal>> upperAgeInput =
            new Input<>("upperAge", "the soft upper bound on the age of the clade", Input.Validate.REQUIRED);
    public Input<RealScalar<NonNegativeReal>> lowerAgeInput =
            new Input<>("lowerAge", "the soft lower bound on the age of the clade", Input.Validate.REQUIRED);
    public Input<RealScalar<UnitInterval>> pCoverageInput =
            new Input<>("confidenceLevel", "the amount of probability mass in the bounds" +
                    "default value (0.9)", new RealScalarParam<>(0.9, UnitInterval.INSTANCE));

    @Override
    public void initAndValidate() {
        if (taxaInput.get().getTaxonSet().isEmpty()) {
            throw new IllegalArgumentException("CalibrationCladePrior " + getID() + " must contain at least one taxon.");
        }
        double p = pCoverageInput.get().get();
        if ((p < 0.0) || (p > 1.0)) {
            throw new IllegalArgumentException("confidenceLevel (" + p + ") should be between 0.0 and 1.0");
        }
        double t_lo = lowerAgeInput.get().get();
        double t_hi = upperAgeInput.get().get();
        if (t_hi < t_lo) {
            throw new IllegalArgumentException("lowerAge (" + t_lo + ") should be less than upperAge (" + t_hi + ")");
        }
        if ((t_lo < 0.0) || (t_hi < 0.0)) {
            throw new IllegalArgumentException("lowerAge should be greater than 0.0: lowerAge=(" + t_lo + "), upperAge=(" + t_hi + ")");
        }
    }

    public TaxonSet getTaxa() {
        return taxaInput.get();
    }

    public boolean isOverlapEdge; // true if overlaps parent interval

    // Fitted parameters
    public double mu;        // log-mean (lognormal)
    public double sigma2;    // log-variance (lognormal)

    double mEdge;
    double vEdge;            // edge log-mean and log-variance increments

    public double alpha;     // Beta alpha
    public double beta;      // Beta beta

    public double getUpper() {
        return upperAgeInput.get().get();
    }

    public double getLower() {
        return lowerAgeInput.get().get();
    }

    public double getCoverage() {
        return pCoverageInput.get().get();
    }
}
