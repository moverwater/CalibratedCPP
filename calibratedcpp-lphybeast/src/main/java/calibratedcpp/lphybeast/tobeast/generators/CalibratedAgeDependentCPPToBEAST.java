package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.type.RealScalar;
import calibratedcpp.CalibratedAgeDependentBirthDeathModel;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.tree.CalibratedAgeDependentCPPTree;
import calibration.CalibrationClade;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedAgeDependentCPPToBEAST
        implements GeneratorToBEAST<CalibratedAgeDependentCPPTree, CalibratedAgeDependentBirthDeathModel> {

    @Override
    public CalibratedAgeDependentBirthDeathModel generatorToBEAST(
            CalibratedAgeDependentCPPTree generator, BEASTInterface value, BEASTContext context) {

        List<TaxonSet> taxonSets = new ArrayList<>();
        CalibratedAgeDependentBirthDeathModel model = new CalibratedAgeDependentBirthDeathModel();
        model.setInputValue("tree", value);

        boolean rootConditioned = generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        model.setInputValue("conditionOnCalibrations", true); // users can override in XML

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(
                    generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        model.setInputValue("birthRate", context.getAsRealScalar(generator.getBirthRate()));
        model.setInputValue("rho",       context.getAsRealScalar(generator.getSamplingProb()));

        model.setInputValue("lifetimeDistribution",
                buildLifetimeDistribution(generator.getLifetime(), context));

        // Remove the 'lifetime' state node and its auto-generated prior; the distribution
        // parameters are wired directly into lifetimeDistribution instead.
        Value<Number> lifetimeValue = generator.getLifetime();
        BEASTInterface lifetimeBEAST = context.getBEASTObject(lifetimeValue);
        if (lifetimeBEAST != null) context.removeBEASTObject(lifetimeBEAST);
        if (lifetimeValue.getGenerator() != null) {
            BEASTInterface lifetimePrior = context.getBEASTObject(lifetimeValue.getGenerator());
            if (lifetimePrior != null) context.removeBEASTObject(lifetimePrior);
        }

        // calibration clades
        List<CalibrationClade> calibrations = new ArrayList<>();
        Calibration[] calibrationsFromGenerator = generator.getCalibrations().value().getCalibrationArray();
        for (Calibration calibration : calibrationsFromGenerator) {
            CalibrationClade calibrationClade = new CalibrationClade();
            TaxonSet taxonSet = getTaxonSet((TreeInterface) value, calibration.getTaxa());
            taxonSets.add(taxonSet);
            calibrationClade.setInputValue("taxa", taxonSet);
            calibrationClade.initAndValidate();
            calibrations.add(calibrationClade);
        }
        model.setInputValue("calibrations", calibrations);
        model.initAndValidate();

        // CalibrationPrior — same pattern as CalibratedCPPToBEAST
        CalibrationPrior calibrationPrior = new CalibrationPrior();
        calibrationPrior.setInputValue("tree", value);
        ConditionedMRCAPrior conditionedMRCAPrior =
                (ConditionedMRCAPrior) generator.getCalibrations().getInputs().get(0);
        Value<Double[]> upperBoundsInput = (Value<Double[]>) conditionedMRCAPrior.getParams().get("upperBounds");
        Value<Double[]> lowerBoundsInput = (Value<Double[]>) conditionedMRCAPrior.getParams().get("lowerBounds");
        Value<Double> covInput = conditionedMRCAPrior.getParams().get("p") != null
                ? (Value<Double>) conditionedMRCAPrior.getParams().get("p")
                : new Value<>("", 0.9);
        List<CalibrationCladePrior> cladeInput = buildCalibrationCladePriors(
                covInput, calibrationsFromGenerator, upperBoundsInput, lowerBoundsInput, taxonSets);
        calibrationPrior.setInputValue("calibration", cladeInput);
        calibrationPrior.initAndValidate();
        context.addBEASTObject(calibrationPrior, conditionedMRCAPrior);

        return model;
    }

    /**
     * Inspects the generator of the LPhy {@code lifetime} value and constructs the
     * corresponding BEAST {@link ScalarDistribution} with its parameters wired from
     * the BEAST context. Named LPhy parameters (e.g. {@code shape.lifetime ~ LogNormal(...)})
     * become live BEAST state nodes; anonymous constants are inlined.
     *
     * <p>Supported LPhy distributions:
     * <ul>
     *   <li>{@code Gamma(shape, scale)}      → BEAST {@code Gamma(alpha=shape, theta=scale)}</li>
     *   <li>{@code LogNormal(meanlog, sdlog)} → BEAST {@code LogNormal(M=meanlog, S=sdlog)}</li>
     *   <li>{@code Exp(mean)}                → BEAST {@code Exponential(mean=mean)}</li>
     *   <li>constant {@code v}               → BEAST {@code Uniform(lower=v−0.01, upper=v+0.01)}</li>
     * </ul>
     */
    private ScalarDistribution buildLifetimeDistribution(Value<Number> lifetime, BEASTContext context) {
        Generator<?> gen = lifetime.getGenerator();
        if (gen == null) {
            double v = lifetime.value().doubleValue();
            return buildUniform(v - 0.01, v + 0.01);
        }

        if (gen instanceof lphy.base.distribution.Gamma g) {
            return buildGamma(g, context);
        } else if (gen instanceof lphy.base.distribution.LogNormal ln) {
            return buildLogNormal(ln, context);
        } else if (gen instanceof lphy.base.distribution.Exp e) {
            return buildExponential(e, context);
        } else {
            throw new IllegalArgumentException(
                    "Unsupported lifetime distribution: " + gen.getClass().getSimpleName()
                    + ". Supported distributions: Gamma, LogNormal, Exp.");
        }
    }

    private Uniform buildUniform(double lower, double upper) {
        Uniform u = new Uniform();
        u.setInputValue("lower", new RealScalarParam<>(lower, Real.INSTANCE));
        u.setInputValue("upper", new RealScalarParam<>(upper, Real.INSTANCE));
        u.initAndValidate();
        return u;
    }

    private Gamma buildGamma(lphy.base.distribution.Gamma g, BEASTContext context) {
        Gamma beastGamma = new Gamma();
        beastGamma.setInputValue("alpha", context.getAsRealScalar(g.getShape()));
        beastGamma.setInputValue("theta", context.getAsRealScalar(g.getScale()));
        beastGamma.initAndValidate();
        return beastGamma;
    }

    private LogNormal buildLogNormal(lphy.base.distribution.LogNormal ln, BEASTContext context) {
        LogNormal beastLN = new LogNormal();
        beastLN.setInputValue("M", context.getAsRealScalar(ln.getMeanLog()));
        beastLN.setInputValue("S", context.getAsRealScalar(ln.getSDLog()));
        beastLN.initAndValidate();
        return beastLN;
    }

    private Exponential buildExponential(lphy.base.distribution.Exp e, BEASTContext context) {
        Exponential beastExp = new Exponential();
        beastExp.setInputValue("mean", context.getAsRealScalar(
                (Value<?>) e.getParams().get("mean")));
        beastExp.initAndValidate();
        return beastExp;
    }

    private List<CalibrationCladePrior> buildCalibrationCladePriors(
            Value<Double> covInput, Calibration[] calibrationsFromGenerator,
            Value<Double[]> upperBoundsInput, Value<Double[]> lowerBoundsInput,
            List<TaxonSet> taxonSets) {

        RealScalar<UnitInterval> confidenceLevel =
                new RealScalarParam<>(covInput.value(), UnitInterval.INSTANCE);
        List<CalibrationCladePrior> cladeInput = new ArrayList<>();
        for (int i = 0; i < calibrationsFromGenerator.length; i++) {
            RealScalar<NonNegativeReal> upperAge =
                    new RealScalarParam<>(upperBoundsInput.value()[i], NonNegativeReal.INSTANCE);
            RealScalar<NonNegativeReal> lowerAge =
                    new RealScalarParam<>(lowerBoundsInput.value()[i], NonNegativeReal.INSTANCE);
            CalibrationCladePrior prior = new CalibrationCladePrior();
            prior.setInputValue("upperAge",        upperAge);
            prior.setInputValue("lowerAge",        lowerAge);
            prior.setInputValue("confidenceLevel", confidenceLevel);
            prior.setInputValue("taxa",            taxonSets.get(i));
            prior.initAndValidate();
            cladeInput.add(prior);
        }
        return cladeInput;
    }

    @Override
    public Class<CalibratedAgeDependentCPPTree> getGeneratorClass() {
        return CalibratedAgeDependentCPPTree.class;
    }

    @Override
    public Class<CalibratedAgeDependentBirthDeathModel> getBEASTClass() {
        return CalibratedAgeDependentBirthDeathModel.class;
    }
}
