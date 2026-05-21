package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.CalibratedBirthDeathModel;
import calibratedcpp.CalibratedBirthDeathSkylineModel;
import calibratedcpp.SkylineParameter;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import calibratedcpp.lphy.tree.CPPTree;

public class CPPToBEAST implements GeneratorToBEAST<CPPTree, CalibratedBirthDeathSkylineModel> {
    @Override
    public CalibratedBirthDeathSkylineModel generatorToBEAST(CPPTree generator, BEASTInterface value, BEASTContext context) {
        CalibratedBirthDeathSkylineModel cpp = new CalibratedBirthDeathSkylineModel();
        cpp.setInputValue("tree", value);

        if (generator.getRootAge() != null) {
            cpp.setInputValue("conditionOnRoot", true);
            cpp.setInputValue("origin", context.getAsRealScalar(generator.getRootAge()));
        } else {
            cpp.setInputValue("origin", new RealScalarParam<>(generator.getConditionAge(), PositiveReal.INSTANCE));
            cpp.setInputValue("conditionOnRoot", false);
        }

        if (generator.getBirthRate() != null) {
            SkylineParameter b = new SkylineParameter();
            b.setInputValue("values", context.getAsRealScalar(generator.getBirthRate()));
            b.initAndValidate();
            cpp.setInputValue("birthRate", b);
        }

        if (generator.getDeathRate() != null) {
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealScalar(generator.getDeathRate()));
            cpp.setInputValue("deathRate", d);
        }

        if (generator.getTurnover() != null) {
            SkylineParameter t = new SkylineParameter();
            t.setInputValue("values", context.getAsRealScalar(generator.getTurnover()));
            t.initAndValidate();
            cpp.setInputValue("turnover", t);
        }

        if (generator.getDiversificationRate() != null) {
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealScalar(generator.getDiversificationRate()));
            d.initAndValidate();
            cpp.setInputValue("diversificationRate", d);
        }

        cpp.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProbability()));

        cpp.initAndValidate();
        return cpp;
    }

    @Override
    public Class<CPPTree> getGeneratorClass() {
        return CPPTree.class;
    }

    @Override
    public Class<CalibratedBirthDeathSkylineModel> getBEASTClass() {
        return CalibratedBirthDeathSkylineModel.class;
    }
}
