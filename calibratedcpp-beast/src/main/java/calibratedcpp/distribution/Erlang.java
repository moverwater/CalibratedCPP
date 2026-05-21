package calibratedcpp.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.IntScalar;
import beast.base.spec.type.RealScalar;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.GammaDistribution;

import java.util.List;

@Description("Erlang distribution: Gamma with a positive-integer shape parameter k and scale θ. " +
        "Density: x^(k-1) * exp(-x/θ) / (θ^k * (k-1)!)  for x > 0.")
public class Erlang extends ScalarDistribution<RealScalar<PositiveReal>, Double> {

    public final Input<IntScalar<PositiveInt>> shapeInput = new Input<>("shape",
            "Erlang shape parameter k (positive integer).");
    public final Input<RealScalar<PositiveReal>> scaleInput = new Input<>("scale",
            "Erlang scale parameter θ (mean = k·θ).");

    private GammaDistribution dist = GammaDistribution.of(1.0, 1.0);
    private ContinuousDistribution.Sampler sampler;

    public Erlang() {}

    @Override
    public void initAndValidate() {
        refresh();
        super.initAndValidate();
    }

    @Override
    public void refresh() {
        int    k     = shapeInput.get().get();
        double theta = scaleInput.get().get();
        if (dist.getShape() != k || dist.getScale() != theta) {
            dist    = GammaDistribution.of(k, theta);
            sampler = null;   // invalidate cached sampler
        }
    }

    @Override
    protected GammaDistribution getApacheDistribution() {
        refresh();
        return dist;
    }

    @Override
    public List<Double> sample() {
        if (sampler == null)
            sampler = dist.createSampler(rng);
        return List.of(sampler.sample());
    }
}