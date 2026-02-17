package calibratedcpp.lphy.spi;

import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.tree.CPPTree;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import calibratedcpp.lphy.util.TruncatedLogNormal;
import lphy.base.spi.LPhyBaseImpl;
import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;

import java.util.Arrays;
import java.util.List;

public class CalibratedcppImpl extends LPhyBaseImpl {
    public CalibratedcppImpl() {}

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList(
            CPPTree.class, CalibratedCPPTree.class, TruncatedLogNormal.class, ConditionedMRCAPrior.class
        );
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Arrays.asList();
    }

    @Override
    public String getExtensionName() {
        return "calibratedcpp lphy library";
    }
}
