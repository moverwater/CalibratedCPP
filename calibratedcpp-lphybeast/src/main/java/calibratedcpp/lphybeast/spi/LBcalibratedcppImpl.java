package calibratedcpp.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.util.TruncatedLogNormal;
import calibratedcpp.lphybeast.tobeast.generators.CPPToBEAST;
import calibratedcpp.lphybeast.tobeast.generators.CalibratedCPPToBEAST;
import jebl.evolution.sequences.SequenceType;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;

import java.util.List;
import java.util.Map;

public class LBcalibratedcppImpl implements LPhyBEASTExt {
    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return List.of(
                CalibratedCPPToBEAST.class, CPPToBEAST.class
                //CalibratedCPPToBEAST_broken.class, CPPToBEAST.class
        );
    }

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return List.of();
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return List.of(
                TruncatedLogNormal.class, ConditionedMRCAPrior.class
        );
    }

    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        return Map.of();
    }

    @Override
    public List<Class> getExcludedValueType() {
        return List.of(
                CalibrationArray.class
        );
    }
}
