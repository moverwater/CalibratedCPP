package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Uniform;
import calibratedcpp.lphy.prior.UniformMRCA;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

/**
 * Converts a single {@link UniformMRCA} calibration into a plain, independently-bounded BEAST
 * {@code MRCAPrior(monophyletic=true, distr=Uniform(lower,upper))}.
 *
 * <p>Registered in {@code LBcalibratedcppImpl.getExcludedGenerator()} rather than
 * {@code getGeneratorToBEASTs()} — {@code UniformMRCA} produces a {@code Calibration}, which
 * (like {@code ConditionedMRCAPrior}'s output) has no independent BEAST representation of its
 * own, and its {@code taxonset} must reuse the exact same {@link TaxonSet} instance the tree
 * model's {@code calibrations} list uses for the same clade. So this class is invoked directly
 * by {@code CalibratedCPPToBEAST}/{@code CalibratedAgeDependentCPPToBEAST} rather than through
 * lphybeast's generic auto-traversal, mirroring how {@code ConditionedMRCAPrior} is handled.
 */
public class UniformMRCAToBEAST implements GeneratorToBEAST<UniformMRCA, MRCAPrior> {

    @Override
    public MRCAPrior generatorToBEAST(UniformMRCA generator, BEASTInterface treeValue, BEASTContext context) {
        TaxonSet taxonSet = getTaxonSet((TreeInterface) treeValue, generator.getTaxa().value());
        return generatorToBEAST(generator, treeValue, taxonSet, context);
    }

    /** Variant that reuses an already-built TaxonSet, so callers can share it with other clades. */
    public MRCAPrior generatorToBEAST(UniformMRCA generator, BEASTInterface treeValue, TaxonSet taxonSet, BEASTContext context) {
        Uniform uniform = new Uniform();
        uniform.setInputValue("lower", context.getAsRealScalar(generator.getLower()));
        uniform.setInputValue("upper", context.getAsRealScalar(generator.getUpper()));
        uniform.initAndValidate();

        MRCAPrior mrcaPrior = new MRCAPrior();
        mrcaPrior.setInputValue("tree", treeValue);
        mrcaPrior.setInputValue("taxonset", taxonSet);
        mrcaPrior.setInputValue("monophyletic", true);
        mrcaPrior.setInputValue("distr", uniform);
        mrcaPrior.initAndValidate();

        context.addBEASTObject(mrcaPrior, generator);
        context.addExtraLoggable(mrcaPrior);
        return mrcaPrior;
    }

    @Override
    public Class<UniformMRCA> getGeneratorClass() {
        return UniformMRCA.class;
    }

    @Override
    public Class<MRCAPrior> getBEASTClass() {
        return MRCAPrior.class;
    }
}
