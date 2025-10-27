package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Prior;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.model.BirthDeathModel;
import lphy.base.evolution.birthdeath.CalibratedCPPTree;
import lphy.core.model.BasicFunction;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedCPPToBEAST implements GeneratorToBEAST<CalibratedCPPTree, CalibratedCoalescentPointProcess> {
    @Override
    public CalibratedCoalescentPointProcess generatorToBEAST(CalibratedCPPTree generator, BEASTInterface value, BEASTContext context) {
        CalibratedCoalescentPointProcess calibratedCPP = new CalibratedCoalescentPointProcess();
        calibratedCPP.setInputValue("tree", value);
        boolean rootConditioned = generator.getRootCondition();
        calibratedCPP.setInputValue("conditionOnRoot", rootConditioned);

        if (rootConditioned){
            // if rootAge is given
            if (generator.getRootAge() != null){
                // if root age is passed in, then use it as real param
                calibratedCPP.setInputValue("origin", context.getAsRealParameter(generator.getRootAge()));
            }
            // otherwise calibration is hidden in clade calibrations, done when looping calibrations

        } else {
            if (generator.getStemAge() != null){
                calibratedCPP.setInputValue("origin", context.getAsRealParameter(generator.getStemAge()));
            }
            // else leave it empty
        }

        // get tree model
        BirthDeathModel treeModel = new BirthDeathModel();
        treeModel.setInputValue("birthRate", context.getAsRealParameter(generator.getBirthRate()));
        treeModel.setInputValue("deathRate", context.getAsRealParameter(generator.getDeathRate()));
        treeModel.setInputValue("rho", context.getAsRealParameter(generator.getSamplingProb()));
        treeModel.initAndValidate();

        calibratedCPP.setInputValue("treeModel", treeModel);

        // get clade calibrations
        int n = generator.getN().value();
        List<CalibrationPoint> calibrations = new ArrayList<>();
        String[][] cladeNames = generator.getCladeTaxa().value();
        Value<Number[]> cladeAges = generator.getCladeAge();
        // TODO: only support distribution generated ages, fixed values will throw error
        BasicFunction tmp = (BasicFunction) cladeAges.getInputs().get(0);

        for (int i = 0; i < cladeNames.length; i++) {
            // get age distribution

            Value cladeAgeValue = tmp.getParams().get(String.valueOf(i));
            BEASTInterface beastCladeAgeValue = context.getBEASTObject(cladeAgeValue);
            context.removeBEASTObject(beastCladeAgeValue);

            if (cladeNames[i].length != n) {
                // get taxon set
                String[] cladeName = new String[cladeNames[i].length];
                int index = 0;
                for (String name : cladeNames[i]) {
                    cladeName[index++] = name;
                }

                TaxonSet cladeTaxonSet = getTaxonSet((TreeInterface) value, cladeName);

                // remove age
                Generator cladePriorGenerator = cladeAgeValue.getGenerator();
                context.removeBEASTObject(beastCladeAgeValue);

                Prior cladeCalibrationPrior = (Prior) context.getBEASTObject(cladePriorGenerator);
                ParametricDistribution calibrationDistribution = cladeCalibrationPrior.distInput.get();
                // remove the clade mrca age in the prior section
                context.removeBEASTObject(cladeCalibrationPrior);

                CalibrationPoint calibrationPoint = new CalibrationPoint();
                calibrationPoint.setInputValue("taxonset", cladeTaxonSet);
                calibrationPoint.setInputValue("distr", calibrationDistribution);
                calibrationPoint.initAndValidate();
                calibrations.add(calibrationPoint);
            } else {
                calibratedCPP.setInputValue("origin", beastCladeAgeValue);
            }
        }

        calibratedCPP.setInputValue("calibrations", calibrations);
        // TODO: make conditionOnCalibrations work with LphyBeast command line options
        calibratedCPP.initAndValidate();

        return calibratedCPP;
    }

    @Override
    public Class<CalibratedCPPTree> getGeneratorClass() {
        return CalibratedCPPTree.class;
    }

    @Override
    public Class<CalibratedCoalescentPointProcess> getBEASTClass() {
        return CalibratedCoalescentPointProcess.class;
    }
}
