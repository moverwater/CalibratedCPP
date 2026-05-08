package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import calibratedcpp.CalibratedBirthDeathModel;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibration.CalibrationClade;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedCPPToBEAST implements GeneratorToBEAST<CalibratedCPPTree, CalibratedBirthDeathModel> {
    @Override
    public CalibratedBirthDeathModel generatorToBEAST(CalibratedCPPTree generator, BEASTInterface value, BEASTContext context) {
        List<TaxonSet> taxonSets = new ArrayList<>();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.setInputValue("tree", value);
        boolean rootConditioned = generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        model.setInputValue("conditionOnCalibrations", true); // users should change this in XML

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        if (generator.getBirthRate() != null)
            model.setInputValue("birthRate", context.getAsRealScalar(generator.getBirthRate()));

        if (generator.getDeathRate() != null)
            model.setInputValue("deathRate", context.getAsRealScalar(generator.getDeathRate()));

        if (generator.getTurnover() != null)
            model.setInputValue("turnover", context.getAsRealScalar(generator.getTurnover()));

        if (generator.getDiversificationRate() != null)
            model.setInputValue("diversificationRate", context.getAsRealScalar(generator.getDiversificationRate()));

        model.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProb()));

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
        ConditionedMRCAPrior conditionedMRCAPrior = (ConditionedMRCAPrior) generator.getCalibrations().getInputs().get(0);
        Value<Double[]> upperBoundsInput = (Value<Double[]>) conditionedMRCAPrior.getParams().get("upperBounds");
        Value<Double[]> lowerBoundsInput = (Value<Double[]>) conditionedMRCAPrior.getParams().get("lowerBounds");
        Value<Double> covInput;
        if (conditionedMRCAPrior.getParams().get("p") != null) {
            covInput = (Value<Double>) conditionedMRCAPrior.getParams().get("p");
        } else {
            covInput = new Value<>("", 0.9);
        }
        List<CalibrationCladePrior> cladeInput = getCalibrationCladePriors(covInput, calibrationsFromGenerator, upperBoundsInput, lowerBoundsInput, taxonSets);
        calibrationPrior.setInputValue("calibration", cladeInput);
        calibrationPrior.initAndValidate();

        context.addBEASTObject(calibrationPrior, conditionedMRCAPrior);

        return model;
    }

    private List<CalibrationCladePrior> getCalibrationCladePriors(Value<Double> covInput, Calibration[] calibrationsFromGenerator, Value<Double[]> upperBoundsInput, Value<Double[]> lowerBoundsInput, List<TaxonSet> taxonSets) {
        RealScalar<UnitInterval> confidenceLevel = new RealScalarParam<>(covInput.value(), UnitInterval.INSTANCE);

        List<CalibrationCladePrior> cladeInput = new ArrayList<>();
        for (int i = 0; i < calibrationsFromGenerator.length; i++) {
            Double upper = upperBoundsInput.value()[i];
            Double lower = lowerBoundsInput.value()[i];

            RealScalar<NonNegativeReal> upperAge = new RealScalarParam<>(upper, NonNegativeReal.INSTANCE);
            RealScalar<NonNegativeReal> lowerAge = new RealScalarParam<>(lower, NonNegativeReal.INSTANCE);

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
    public Class<CalibratedBirthDeathModel> getBEASTClass() {
        return CalibratedBirthDeathModel.class;
    }
}
