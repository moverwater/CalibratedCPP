package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.spec.domain.Real;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.OffsetExponentialMRCA;
import calibratedcpp.lphy.prior.UniformMRCA;
import calibratedcpp.lphy.prior.toCalibrationArray;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphy.core.vectorization.array.ArrayFunction;
import lphybeast.BEASTContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Support for calibrations built from independent per-clade {@link UniformMRCA} /
 * {@link OffsetExponentialMRCA} LPhy generators (combined via {@code toArray(calibrations=[cal1,
 * cal2, ...])}, the same array-literal syntax used by {@code ConditionedMRCAPrior(calibrations=
 * [cal1, cal2, ...])}), as opposed to the joint calibration density produced by
 * {@link calibratedcpp.lphy.prior.ConditionedMRCAPrior}. The array can freely mix both generator
 * types -- each element owns its own age-distribution shape (Uniform vs. offset-exponential) and
 * is dispatched to the matching converter independently.
 */
public class MRCAPriorCalibrationUtils {

    /**
     * @return true if {@code calibrationsValue}'s generator is not {@code ConditionedMRCAPrior}
     *         (i.e. it was built via {@code toArray(calibrations=[...])} over independent
     *         {@code UniformMRCA}/{@code OffsetExponentialMRCA} calls instead).
     */
    public static boolean isIndependentMRCAPriorSource(Value<CalibrationArray> calibrationsValue) {
        return calibrationsValue.getGenerator() instanceof toCalibrationArray;
    }

    /**
     * Walks the {@code toArray(calibrations=[...])} array literal feeding {@code calibrationsValue}
     * and returns the individual per-clade calibration generators (each a {@code UniformMRCA} or
     * {@code OffsetExponentialMRCA}) in the same order their calibrations appear in the resulting
     * {@link CalibrationArray} — so index i here lines up with index i of
     * {@code calibrationsValue.value().getCalibrationArray()} and with any TaxonSet list built
     * from that same array.
     */
    public static List<Generator<?>> collectIndependentCalibrationGenerators(Value<?> calibrationsValue) {
        List<Generator<?>> out = new ArrayList<>();
        Generator<?> gen = calibrationsValue.getGenerator();
        if (!(gen instanceof toCalibrationArray tca)) {
            throw new IllegalArgumentException(
                    "Expected calibrations to come from toArray(calibrations=[...]), got: "
                            + (gen == null ? "constant value" : gen.getClass().getSimpleName()));
        }
        Value<?> arrayValue = tca.getParams().get(toCalibrationArray.calibrationsParamName);
        Generator<?> arrayGen = arrayValue.getGenerator();
        if (!(arrayGen instanceof ArrayFunction<?> arrayFunction)) {
            throw new IllegalArgumentException(
                    "Expected an array literal of calibrations (e.g. [cal1, cal2, ...]), got: "
                            + (arrayGen == null ? "constant value" : arrayGen.getClass().getSimpleName()));
        }
        for (Value<?> element : arrayFunction.getValues()) {
            Generator<?> elementGen = element.getGenerator();
            if (!(elementGen instanceof UniformMRCA) && !(elementGen instanceof OffsetExponentialMRCA)) {
                throw new IllegalArgumentException(
                        "Expected each calibration in the array to come from UniformMRCA or OffsetExponentialMRCA, got: "
                                + (elementGen == null ? "constant value" : elementGen.getClass().getSimpleName()));
            }
            out.add(elementGen);
        }
        return out;
    }

    /**
     * Builds one plain BEAST {@code MRCAPrior} per calibration generator (dispatched to
     * {@code UniformMRCAToBEAST} or {@code OffsetExponentialMRCAToBEAST} by its actual type),
     * reusing the caller's already-built {@link TaxonSet} for each clade (index-aligned with
     * {@code calibrationGenerators}) so the {@code taxonset} reference matches the one used in
     * the tree model's {@code calibrations} list, rather than building a second, duplicate
     * TaxonSet for the same clade.
     */
    public static void buildIndependentMRCAPriors(
            List<Generator<?>> calibrationGenerators, List<TaxonSet> taxonSets,
            BEASTInterface treeValue, BEASTContext context) {

        UniformMRCAToBEAST uniformConverter = new UniformMRCAToBEAST();
        OffsetExponentialMRCAToBEAST offsetExponentialConverter = new OffsetExponentialMRCAToBEAST();
        for (int i = 0; i < calibrationGenerators.size(); i++) {
            Generator<?> gen = calibrationGenerators.get(i);
            if (gen instanceof UniformMRCA uniformMRCA) {
                uniformConverter.generatorToBEAST(uniformMRCA, treeValue, taxonSets.get(i), context);
            } else if (gen instanceof OffsetExponentialMRCA offsetExponentialMRCA) {
                offsetExponentialConverter.generatorToBEAST(offsetExponentialMRCA, treeValue, taxonSets.get(i), context);
            } else {
                throw new IllegalArgumentException("Unsupported calibration generator: " + gen.getClass().getSimpleName());
            }
        }
    }

    /**
     * Builds a plain BEAST {@code MRCAPrior(monophyletic=true, distr=Uniform(lower,upper))} from
     * raw double bounds, for calibration sources (e.g. {@code ConditionedMRCAPrior}) that only
     * expose bounds as plain numbers rather than as LPhy {@code Value}s.
     */
    public static MRCAPrior buildBoundedMRCAPrior(BEASTInterface treeValue, TaxonSet taxonSet, double lower, double upper) {
        Uniform uniform = new Uniform();
        uniform.setInputValue("lower", new RealScalarParam<>(lower, Real.INSTANCE));
        uniform.setInputValue("upper", new RealScalarParam<>(upper, Real.INSTANCE));
        uniform.initAndValidate();

        MRCAPrior mrcaPrior = new MRCAPrior();
        mrcaPrior.setInputValue("tree", treeValue);
        mrcaPrior.setInputValue("taxonset", taxonSet);
        mrcaPrior.setInputValue("monophyletic", true);
        mrcaPrior.setInputValue("distr", uniform);
        mrcaPrior.initAndValidate();
        return mrcaPrior;
    }
}
