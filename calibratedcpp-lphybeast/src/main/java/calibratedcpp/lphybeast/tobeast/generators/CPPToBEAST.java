package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.CalibratedBirthDeathModel;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import calibratedcpp.lphy.tree.CPPTree;

public class CPPToBEAST implements GeneratorToBEAST<CPPTree, CalibratedBirthDeathModel> {
    @Override
    public CalibratedBirthDeathModel generatorToBEAST(CPPTree generator, BEASTInterface value, BEASTContext context) {
        CalibratedBirthDeathModel cpp = new CalibratedBirthDeathModel();
        cpp.setInputValue("tree", value);

        if (generator.getRootAge() != null) {
            cpp.setInputValue("conditionOnRoot", true);
            cpp.setInputValue("origin", context.getAsRealScalar(generator.getRootAge()));
        } else {
            cpp.setInputValue("origin", new RealScalarParam<>(generator.getConditionAge(), PositiveReal.INSTANCE));
            cpp.setInputValue("conditionOnRoot", false);
        }

        if (generator.getBirthRate() != null)
            cpp.setInputValue("birthRate", context.getAsRealScalar(generator.getBirthRate()));

        if (generator.getDeathRate() != null)
            cpp.setInputValue("deathRate", context.getAsRealScalar(generator.getDeathRate()));

        if (generator.getTurnover() != null)
            cpp.setInputValue("turnover", context.getAsRealScalar(generator.getTurnover()));

        if (generator.getDiversificationRate() != null)
            cpp.setInputValue("diversificationRate", context.getAsRealScalar(generator.getDiversificationRate()));

        cpp.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProbability()));

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
