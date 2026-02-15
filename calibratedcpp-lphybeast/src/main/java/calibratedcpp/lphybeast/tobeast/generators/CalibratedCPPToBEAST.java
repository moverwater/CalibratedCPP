package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.model.BirthDeathModel;
import calibration.CalibrationClade;
import lphy.base.evolution.birthdeath.CalibratedCPPTree;
import lphy.core.model.BasicFunction;
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
        List<CalibrationClade> calibrations = new ArrayList<>();
        String[][] cladeNames = generator.getCladeTaxa().value();
        Value<Number[]> cladeAges = generator.getCladeAge();

        // check if age is generated from distributions
        if (cladeAges.getInputs().size() != 0 ) {
            BasicFunction tmp = (BasicFunction) cladeAges.getInputs().get(0);
            for (int i = 0; i < cladeNames.length; i++) {
                if (cladeNames[i].length != n) {
                    // get taxon set
                    String[] cladeName = new String[cladeNames[i].length];
                    int index = 0;
                    for (String name : cladeNames[i]) {
                        cladeName[index++] = name;
                    }

                    TaxonSet cladeTaxonSet = getTaxonSet((TreeInterface) value, cladeName);

                    CalibrationClade clade = new CalibrationClade();
                    clade.setInputValue("taxa", cladeTaxonSet);

                    RealParameter age = context.getAsRealParameter(tmp.getParams().get(String.valueOf(i)));

                    clade.setInputValue("tmrca", age);

                    clade.initAndValidate();
                    calibrations.add(clade);

                } else {
                    RealParameter age = context.getAsRealParameter(tmp.getParams().get(String.valueOf(i)));
                    calibratedCPP.setInputValue("origin", age);
                }

                // remove the separate beast object of clade ages
                Value cladeAgeValue = tmp.getParams().get(String.valueOf(i));
                BEASTInterface beastCladeAgeValue = context.getBEASTObject(cladeAgeValue);
                context.removeBEASTObject(beastCladeAgeValue);
            }
        } else {
            // fixed values as input
            for (int i = 0; i < cladeNames.length; i++) {
                if (cladeNames[i].length != n) {
                    // get taxon set
                    String[] cladeName = new String[cladeNames[i].length];
                    int index = 0;
                    for (String name : cladeNames[i]) {
                        cladeName[index++] = name;
                    }

                    TaxonSet cladeTaxonSet = getTaxonSet((TreeInterface) value, cladeName);

                    CalibrationClade clade = new CalibrationClade();
                    clade.setInputValue("taxa", cladeTaxonSet);

                    RealParameter age = new RealParameter(new Double[]{cladeAges.value()[i].doubleValue()});
                    age.initAndValidate();
                    clade.setInputValue("tmrca", age);

                    clade.initAndValidate();
                    calibrations.add(clade);
                } else {
                    RealParameter age = new RealParameter(new Double[]{cladeAges.value()[i].doubleValue()});
                    age.initAndValidate();
                    calibratedCPP.setInputValue("origin", age);
                }
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
