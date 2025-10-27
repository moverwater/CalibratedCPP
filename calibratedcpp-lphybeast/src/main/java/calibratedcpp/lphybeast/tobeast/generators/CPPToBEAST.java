package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.model.BirthDeathModel;
import lphy.base.evolution.birthdeath.CPPTree;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

public class CPPToBEAST implements GeneratorToBEAST<CPPTree, CalibratedCoalescentPointProcess> {
    @Override
    public CalibratedCoalescentPointProcess generatorToBEAST(CPPTree generator, BEASTInterface value, BEASTContext context) {
        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.setInputValue("tree", value);

        // set origin and condition on root
        // TODO: make it reasonable to let stem age random
        if (generator.getRootAge() != null) {
            cpp.setInputValue("conditionOnRoot", true);
            if (generator.getRootAge().value().doubleValue() == generator.getConditionAge()) {
                cpp.setInputValue("origin", context.getAsRealParameter(generator.getRootAge()));
            } else {
                cpp.setInputValue("origin", generator.getConditionAge());
            }
        } else {
            cpp.setInputValue("origin", generator.getConditionAge());
            cpp.setInputValue("conditionOnRoot", false);
        }

        // get tree model
        BirthDeathModel treeModel = new BirthDeathModel();
        treeModel.setInputValue("birthRate", context.getAsRealParameter(generator.getBirthRate()));
        treeModel.setInputValue("deathRate", context.getAsRealParameter(generator.getDeathRate()));
        treeModel.setInputValue("rho", context.getAsRealParameter(generator.getSamplingProbability()));
        treeModel.initAndValidate();

        cpp.setInputValue("treeModel", treeModel);

        cpp.initAndValidate();
        return cpp;
    }

    @Override
    public Class<CPPTree> getGeneratorClass() {
        return CPPTree.class;
    }

    @Override
    public Class<CalibratedCoalescentPointProcess> getBEASTClass() {
        return CalibratedCoalescentPointProcess.class;
    }
}
