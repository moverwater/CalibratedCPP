package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import calibratedcpp.BirthDeathModel;
import calibratedcpp.lphy.prior.Calibration;
import calibration.CalibrationClade;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import lphy.core.model.BasicFunction;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedCPPToBEAST implements GeneratorToBEAST<CalibratedCPPTree, BirthDeathModel> {
    @Override
    public BirthDeathModel generatorToBEAST(CalibratedCPPTree generator, BEASTInterface value, BEASTContext context) {
        BirthDeathModel model = new BirthDeathModel();
        model.setInputValue("tree", value);
        boolean rootConditioned = generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);

        // if it is not root conditioned, then set the stem age
        // if there is root calibration, then it will be set in initAndValidate
        if (! rootConditioned){
           model.setInputValue("origin", context.getAsRealParameter(generator.getOrigin()));
        }

        // get tree model
        BirthDeathModel treeModel = new BirthDeathModel();
        treeModel.setInputValue("birthRate", context.getAsRealParameter(generator.getBirthRate()));
        treeModel.setInputValue("deathRate", context.getAsRealParameter(generator.getDeathRate()));
        treeModel.setInputValue("rho", context.getAsRealParameter(generator.getSamplingProb()));
        treeModel.initAndValidate();

        model.setInputValue("treeModel", treeModel);

        // get clade calibrations
        List<CalibrationClade> calibrations = new ArrayList<>();
        Calibration[] calibrationsFromGenerator = generator.getCalibrations().value();
        for (Calibration calibration : calibrationsFromGenerator) {
            CalibrationClade calibrationClade = new CalibrationClade();
            TaxonSet taxonSet = getTaxonSet((TreeInterface) value, calibration.getTaxa());
            calibrationClade.setInputValue("taxa", taxonSet);
            calibrationClade.initAndValidate();
            calibrations.add(calibrationClade);
        }

        model.setInputValue("calibrations", calibrations);
        model.initAndValidate();

        /*
            map the conditioned mrca prior section
         */
        CalibrationPrior calibrationPrior = new CalibrationPrior();
        calibrationPrior.setInputValue("tree", value);
        // get the calibrations
        List<CalibrationCladePrior> bounds = new ArrayList<>();
        BasicFunction tmp = (BasicFunction) generator.getCalibrations().getInputs().get(0);
        for (int i = 0; i < calibrationsFromGenerator.length; i++) {

        }


        CalibrationCladePrior calibrationCladePrior = new CalibrationCladePrior();


        return model;
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
