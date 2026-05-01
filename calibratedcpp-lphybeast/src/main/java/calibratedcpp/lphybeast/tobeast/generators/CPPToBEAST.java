package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
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

        // set origin and condition on root
        if (generator.getRootAge() != null) {
            cpp.setInputValue("conditionOnRoot", true);
            cpp.setInputValue("origin", context.getAsRealParameter(generator.getRootAge()));
        } else {
            cpp.setInputValue("origin", generator.getConditionAge());
            cpp.setInputValue("conditionOnRoot", false);
        }

        if (generator.getBirthRate() != null){
            SkylineParameter b = new SkylineParameter();
            b.setInputValue("values", context.getAsRealParameter(generator.getBirthRate()));
            b.initAndValidate();
            cpp.setInputValue("birthRate", b);
        }

        if (generator.getDeathRate() != null){
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealParameter(generator.getDeathRate()));
            cpp.setInputValue("deathRate", d);
        }

        if (generator.getTurnover() != null){
            SkylineParameter t = new SkylineParameter();
            t.setInputValue("values", context.getAsRealParameter(generator.getTurnover()));
            t.initAndValidate();
            cpp.setInputValue("turnover", t);
        }

        if (generator.getDiversificationRate() != null){
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealParameter(generator.getDiversificationRate()));
            d.initAndValidate();
            cpp.setInputValue("diversificationRate", d);
        }

        cpp.setInputValue("rho", context.getAsRealParameter(generator.getSamplingProbability()));

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
