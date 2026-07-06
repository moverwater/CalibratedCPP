package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import calibratedcpp.CalibratedBirthDeathSkylineModel;
import calibratedcpp.SkylineParameter;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedCPPToBEAST implements GeneratorToBEAST<CalibratedCPPTree, CalibratedBirthDeathSkylineModel> {
    @Override
    public CalibratedBirthDeathSkylineModel generatorToBEAST(CalibratedCPPTree generator, BEASTInterface value, BEASTContext context) {
        List<TaxonSet> taxonSets = new ArrayList<>();
        CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
        model.setInputValue("tree", value);

        // When no calibrations are provided rootAge drives conditioning, so conditionOnRoot=true.
        boolean hasCalibrations = generator.getCalibrations() != null;
        boolean rootConditioned = !hasCalibrations || generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        model.setInputValue("conditionOnCalibrations", hasCalibrations);

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        if (generator.getBirthRate() != null) {
            SkylineParameter b = new SkylineParameter();
            b.setInputValue("values", context.getAsRealTensor(generator.getBirthRate()));
            b.initAndValidate();
            model.setInputValue("birthRate", b);
        }

        if (generator.getDeathRate() != null) {
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealTensor(generator.getDeathRate()));
            d.initAndValidate();
            model.setInputValue("deathRate", d);
        }

        if (generator.getTurnover() != null) {
            SkylineParameter t = new SkylineParameter();
            t.setInputValue("values", context.getAsRealTensor(generator.getTurnover()));
            t.initAndValidate();
            model.setInputValue("turnover", t);
        }

        if (generator.getDiversificationRate() != null) {
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealTensor(generator.getDiversificationRate()));
            d.initAndValidate();
            model.setInputValue("diversificationRate", d);
        }

        model.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProb()));

        if (!hasCalibrations) {
            // rootAge-only mode: no calibration clades or CalibrationPrior needed.
            model.setInputValue("calibrations", new ArrayList<>());
            model.initAndValidate();
            return model;
        }

        // get clade calibrations
        Calibration[] calibrationsFromGenerator = generator.getCalibrations().value().getCalibrationArray();

        for (Calibration calibration : calibrationsFromGenerator) {
            TaxonSet taxonSet = getTaxonSet((TreeInterface) value, calibration.getTaxa());
            taxonSets.add(taxonSet);
        }

        model.setInputValue("calibrations", taxonSets);
        model.initAndValidate();

        Value<CalibrationArray> calibrationsValue = generator.getCalibrations();
        if (MRCAPriorCalibrationUtils.isIndependentMRCAPriorSource(calibrationsValue)) {
            // calibrations came from array of independent UniformMRCA(taxa=,upper=,lower=)
            // not ConditionedMRCAPrior's joint density
            // build one plain, independently-bounded MRCAPrior per clade instead of a CalibrationPrior.
            MRCAPriorCalibrationUtils.buildIndependentMRCAPriors(
                    MRCAPriorCalibrationUtils.collectUniformMRCAs(calibrationsValue),
                    taxonSets, value, context);
            return model;
        }

        /*
            map the conditioned MRCA prior objects
         */
        CalibrationPrior calibrationPrior = new CalibrationPrior();
        calibrationPrior.setInputValue("tree", value);
        ConditionedMRCAPrior conditionedMRCAPrior = (ConditionedMRCAPrior) calibrationsValue.getInputs().get(0);
        Value<Calibration[]> calibrationSpecsInput = conditionedMRCAPrior.getCalibrations();
        Value<Double> covInput = conditionedMRCAPrior.getCoverage() != null
                ? new Value<>("", conditionedMRCAPrior.getCoverage().value().doubleValue())
                : new Value<>("", 0.9);
        List<CalibrationCladePrior> cladeInput = getCalibrationCladePriors(covInput, calibrationsFromGenerator, calibrationSpecsInput.value(), taxonSets);
        calibrationPrior.setInputValue("calibration", cladeInput);
        calibrationPrior.initAndValidate();

        context.addExtraLoggable(calibrationPrior);
        context.addBEASTObject(calibrationPrior, conditionedMRCAPrior);

        return model;
    }

    private List<CalibrationCladePrior> getCalibrationCladePriors(Value<Double> covInput, Calibration[] calibrationsFromGenerator, Calibration[] calibrationSpecs, List<TaxonSet> taxonSets) {
        RealScalar<UnitInterval> confidenceLevel = new RealScalarParam<>(covInput.value(), UnitInterval.INSTANCE);

        List<CalibrationCladePrior> cladeInput = new ArrayList<>();
        for (int i = 0; i < calibrationsFromGenerator.length; i++) {
            RealScalar<NonNegativeReal> upperAge = new RealScalarParam<>(calibrationSpecs[i].getUpper(), NonNegativeReal.INSTANCE);
            RealScalar<NonNegativeReal> lowerAge = new RealScalarParam<>(calibrationSpecs[i].getLower(), NonNegativeReal.INSTANCE);

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
    public Class<CalibratedBirthDeathSkylineModel> getBEASTClass() {
        return CalibratedBirthDeathSkylineModel.class;
    }
}
