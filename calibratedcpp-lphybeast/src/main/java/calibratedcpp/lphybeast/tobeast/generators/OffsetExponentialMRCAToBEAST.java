package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.OffsetReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.lphy.prior.OffsetExponentialMRCA;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

/**
 * Converts a single {@link OffsetExponentialMRCA} calibration into a plain, independently-bounded
 * BEAST {@code MRCAPrior(monophyletic=true, distr=OffsetReal(offset=offset, distribution=Exponential(mean=mean)))}.
 *
 * <p>Registered in {@code LBcalibratedcppImpl.getExcludedGenerator()} rather than
 * {@code getGeneratorToBEASTs()}, mirroring {@code UniformMRCAToBEAST} exactly -- see that class's
 * javadoc for why (its output has no independent BEAST representation on its own, and its
 * {@code taxonset} must reuse the same {@link TaxonSet} instance the tree model's calibrations
 * list uses for the same clade).
 */
public class OffsetExponentialMRCAToBEAST implements GeneratorToBEAST<OffsetExponentialMRCA, MRCAPrior> {

    @Override
    public MRCAPrior generatorToBEAST(OffsetExponentialMRCA generator, BEASTInterface treeValue, BEASTContext context) {
        TaxonSet taxonSet = getTaxonSet((TreeInterface) treeValue, generator.getTaxa().value());
        return generatorToBEAST(generator, treeValue, taxonSet, context);
    }

    /** Variant that reuses an already-built TaxonSet, so callers can share it with other clades. */
    public MRCAPrior generatorToBEAST(OffsetExponentialMRCA generator, BEASTInterface treeValue, TaxonSet taxonSet, BEASTContext context) {
        // Exponential.mean requires the PositiveReal domain specifically (not the unconstrained
        // Real domain context.getAsRealScalar produces), so it's built explicitly here rather
        // than via the generic Value converter -- same reasoning as MRCAPriorCalibrationUtils'
        // Uniform builder, which also constructs its RealScalarParams directly.
        Exponential exponential = new Exponential();
        exponential.setInputValue("mean", new RealScalarParam<>(generator.getMean().value().doubleValue(), PositiveReal.INSTANCE));
        exponential.initAndValidate();

        OffsetReal offsetReal = new OffsetReal();
        offsetReal.setInputValue("offset", context.getAsRealScalar(generator.getOffset()));
        offsetReal.setInputValue("distribution", exponential);
        offsetReal.initAndValidate();

        MRCAPrior mrcaPrior = new MRCAPrior();
        mrcaPrior.setInputValue("tree", treeValue);
        mrcaPrior.setInputValue("taxonset", taxonSet);
        mrcaPrior.setInputValue("monophyletic", true);
        mrcaPrior.setInputValue("distr", offsetReal);
        mrcaPrior.initAndValidate();

        context.addBEASTObject(mrcaPrior, generator);
        context.addExtraLoggable(mrcaPrior);
        return mrcaPrior;
    }

    @Override
    public Class<OffsetExponentialMRCA> getGeneratorClass() {
        return OffsetExponentialMRCA.class;
    }

    @Override
    public Class<MRCAPrior> getBEASTClass() {
        return MRCAPrior.class;
    }
}
