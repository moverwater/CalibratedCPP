package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import calibratedcpp.CalibratedBirthDeathModel;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import calibratedcpp.lphy.tree.CPPTree;

public class CPPToBEAST implements GeneratorToBEAST<CPPTree, CalibratedBirthDeathModel> {
    @Override
    public CalibratedBirthDeathModel generatorToBEAST(CPPTree generator, BEASTInterface value, BEASTContext context) {
        CalibratedBirthDeathModel cpp = new CalibratedBirthDeathModel();
        cpp.setInputValue("tree", value);

        // set origin and condition on root
        if (generator.getRootAge() != null) {
            cpp.setInputValue("conditionOnRoot", true);
            cpp.setInputValue("origin", context.getAsRealParameter(generator.getRootAge()));
        } else {
            cpp.setInputValue("origin", generator.getConditionAge());
            cpp.setInputValue("conditionOnRoot", false);
        }

        // get tree model
        cpp.setInputValue("birthRate", context.getAsRealParameter(generator.getBirthRate()));
        cpp.setInputValue("deathRate", context.getAsRealParameter(generator.getDeathRate()));
        cpp.setInputValue("rho", context.getAsRealParameter(generator.getSamplingProbability()));

        cpp.initAndValidate();
        return cpp;
    }

    @Override
    public Class<CPPTree> getGeneratorClass() {
        return CPPTree.class;
    }

    @Override
    public Class<CalibratedBirthDeathModel> getBEASTClass() {
        return CalibratedBirthDeathModel.class;
    }
}
