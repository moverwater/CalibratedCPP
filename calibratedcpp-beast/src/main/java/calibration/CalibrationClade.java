package calibration;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.parameter.RealParameter;


/**
 * @author Marcus Overwater
 */

@Description("A calibration clade is a set of taxa and an age for the MRCA of the taxa.")
public class CalibrationClade extends BEASTObject {
    public Input<TaxonSet> taxaInput =
            new Input<>("taxa", "The set of taxa in the clade", Input.Validate.REQUIRED);
    public Input<RealParameter> ageInput =
            new Input<>("tmrca", "The time of the most recent common ancestor in the clade", Input.Validate.OPTIONAL);

    public boolean providedAge;

    @Override
    public void initAndValidate() {
        if (taxaInput.get().getTaxonSet().isEmpty()) {
            throw new IllegalArgumentException("Calibration clade " + getID() + " must contain at least one taxon.");
        }
        if (ageInput.get() != null) {
            if (ageInput.get().getValue() <= 0) {
                throw new IllegalArgumentException("Calibration age must be positive: " + getID());
            }
            providedAge = true;
        } else {
            providedAge = false;
        }
    }

    public TaxonSet getTaxa() {
        return taxaInput.get();
    }

    public RealParameter getAge() {
        return ageInput.get();
    }
}
