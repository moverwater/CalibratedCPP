package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.BirthDeathModel;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibration.CalibrationClade;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedCPPToBEAST implements GeneratorToBEAST<CalibratedCPPTree, BirthDeathModel> {
    List<TaxonSet> taxonSets = new ArrayList<>();
    @Override
    public BirthDeathModel generatorToBEAST(CalibratedCPPTree generator, BEASTInterface value, BEASTContext context) {
        BirthDeathModel model = new BirthDeathModel();
        model.setInputValue("tree", value);
        boolean rootConditioned = generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);

        // if it is not root conditioned, then set the stem age
        // if there is root calibration, then it will be set in initAndValidate
        if (! rootConditioned){
            RealParameter realParameter = new RealParameter(new Double[]{generator.getOrigin().value()});
           model.setInputValue("origin", realParameter);
        }

        // get tree model
        model.setInputValue("birthRate", context.getAsRealParameter(generator.getBirthRate()));
        model.setInputValue("deathRate", context.getAsRealParameter(generator.getDeathRate()));
        model.setInputValue("rho", context.getAsRealParameter(generator.getSamplingProb()));

        // get clade calibrations
        List<CalibrationClade> calibrations = new ArrayList<>();
        Calibration[] calibrationsFromGenerator = generator.getCalibrations().value().getCalibrationArray();

        for (Calibration calibration : calibrationsFromGenerator) {
            CalibrationClade calibrationClade = new CalibrationClade();
            TaxonSet taxonSet = getTaxonSet((TreeInterface) value, calibration.getTaxa());
            taxonSets.add(taxonSet);
            calibrationClade.setInputValue("taxa", taxonSet);
            calibrationClade.initAndValidate();
            calibrations.add(calibrationClade);
        }

        model.setInputValue("calibrations", calibrations);
        model.initAndValidate();

        /*
            map the conditioned MRCA prior objects
         */
        CalibrationPrior calibrationPrior = new CalibrationPrior();
        calibrationPrior.setInputValue("tree", value);
        // get the calibrations
        List<CalibrationCladePrior> bounds = new ArrayList<>();
        // extract bounds and coverage info from generator
        ConditionedMRCAPrior conditionedMRCAPrior = (ConditionedMRCAPrior) generator.getCalibrations().getInputs().get(0);
        Value<Double[]> upperBoundsInput = ( Value<Double[]>) conditionedMRCAPrior.getParams().get("upperBounds");
        Value<Double[]> lowerBoundsInput = ( Value<Double[]>) conditionedMRCAPrior.getParams().get("lowerBounds");
        Value<Double> covInput;
        if (conditionedMRCAPrior.getParams().get("p") != null){
            covInput = ( Value<Double>) conditionedMRCAPrior.getParams().get("p");
        } else {
            covInput = new Value<>("", 0.9);
        }
        List<CalibrationCladePrior> cladeInput = getCalibrationCladePriors(covInput, calibrationsFromGenerator, upperBoundsInput, lowerBoundsInput);
        calibrationPrior.setInputValue("calibration", cladeInput);
        calibrationPrior.initAndValidate();

        // add beast beast object for calibrationPrior
       context.addBEASTObject(calibrationPrior,conditionedMRCAPrior);

        return model;
    }

    private List<CalibrationCladePrior> getCalibrationCladePriors(Value<Double> covInput, Calibration[] calibrationsFromGenerator, Value<Double[]> upperBoundsInput, Value<Double[]> lowerBoundsInput) {
        RealParameter confidenceLevel = new RealParameter(new Double[]{covInput.value()});

        List<CalibrationCladePrior> cladeInput = new ArrayList<>();
        for (int i = 0; i < calibrationsFromGenerator.length; i++) {
            Double upper = upperBoundsInput.value()[i];
            Double lower = lowerBoundsInput.value()[i];

            RealParameter upperAge = new RealParameter(new Double[]{upper});
            RealParameter lowerAge = new RealParameter(new Double[]{lower});

            // construct calibration clade prior
            CalibrationCladePrior calibrationCladePrior = new CalibrationCladePrior();
            calibrationCladePrior.setInputValue("upperAge", upperAge);
            calibrationCladePrior.setInputValue("lowerAge", lowerAge);
            calibrationCladePrior.setInputValue("confidenceLevel", confidenceLevel);
            calibrationCladePrior.setInputValue("taxa", taxonSets.get(i));
            calibrationCladePrior.initAndValidate();
            cladeInput.add(calibrationCladePrior);
        }
        return cladeInput;
    }

    @Override
    public Class<CalibratedCPPTree> getGeneratorClass() {
        return CalibratedCPPTree.class;
    }

    @Override
    public Class<BirthDeathModel> getBEASTClass() {
        return BirthDeathModel.class;
    }
}
