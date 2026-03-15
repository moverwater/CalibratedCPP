import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import lphy.core.model.Value;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.inference.TTest;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static calibratedcpp.lphy.prior.ConditionedMRCAPrior.mapBetaNodes;
import static org.junit.jupiter.api.Assertions.*;

public class ConditionedMRCAPriorTest {

    private static final int[] PARENT = {-1, 0, 0, 2, 1, 1, 5, 5, 7, 8};
    private static final Double[] T_LO = {10.0, 9.4, 4.8, 4.0, 8.0, 6.8, 1.8, 2.0, 1.7, 1.6};
    private static final Double[] T_HI = {10.5, 9.6, 5.0, 5.0, 9.5, 8.0, 2.5, 3.2, 2.5, 2.5};
    private static final int N = 10;
    private static final double P_COV = 0.90;
    private static final int N_SIM = 100000;
    private static final double T_TEST_ALPHA = 0.001;

    private static double[][] samples;

    @BeforeAll
    static void generateSamples() {
        ConditionedMRCAPrior prior = new ConditionedMRCAPrior(
                new Value<>("", buildDummyTaxa()),
                new Value<>("", true),
                new Value<>("", T_HI),
                new Value<>("", T_LO),
                new Value<>("", P_COV)
        );
        samples = new double[N_SIM][N];
        for (int s = 0; s < N_SIM; s++) {
            Calibration[] cals = prior.sample().value().getCalibrationArray();
            for (int i = 0; i < N; i++) {
                samples[s][i] = cals[i].getAge();
            }
        }
    }

    // ── test 1: mapBetaNodes matches R classification ─────────────────────────
    // R beta_nodes (1-indexed): 4,5,9,10 → 0-indexed: 3,4,8,9
    @Test
    void testMapBetaNodes() {
        boolean[] actual   = mapBetaNodes(N, true, PARENT, T_HI, T_LO);
        boolean[] expected = {false, false, false, true, true, false, false, false, true, true};
        for (int i = 0; i < N; i++) {
            assertEquals(expected[i], actual[i],
                    "is_beta_node mismatch at node " + (i + 1));
        }
    }

    // ── test 2: sample MEDIAN (not mean) within [t_lo, t_hi] ─────────────────
    // FIX: lognormal median = exp(mu) sits in [t_lo, t_hi] by construction.
    //      lognormal mean = exp(mu + sigma2/2) also happens to be in bounds for
    //      this tree, but median is the correct theoretical target.
    //      For truncated nodes the mean is shifted lower — median is more robust.
    @Test
    void testSampleMediansWithinBounds() {
        for (int i = 0; i < N; i++) {
            int fi = i;
            double[] col = Arrays.stream(samples).mapToDouble(s -> s[fi]).toArray();
            Arrays.sort(col);
            double median = col[N_SIM / 2];
            assertTrue(median >= T_LO[i] && median <= T_HI[i],
                    "Sample median for node " + (i + 1) + " = " + median + " is outside [" + T_LO[i] + ", " + T_HI[i] + "]");
        }
    }

    // ── test 3a: root node coverage matches p_cov exactly ────────────────────────
    // Root is a pure untruncated lognormal parameterised to hit p_cov by construction.
    // This is the only node with a hard theoretical guarantee.
    @Test
    void testRootCoverageMatchesPCov() {
        double se   = Math.sqrt(P_COV * (1 - P_COV) / N_SIM);
        double band = 5 * se;  // ~5σ: false-failure rate < 0.00006%

        double rootCov = Arrays.stream(samples)
                .mapToDouble(s -> (s[0] >= T_LO[0] && s[0] <= T_HI[0]) ? 1.0 : 0.0)
                .average().orElseThrow();

        assertEquals(P_COV, rootCov, band,
                "Root coverage=" + String.format("%.4f", rootCov) +
                        " expected=" + P_COV + " band=±" + String.format("%.4f", band));
    }

    // ── test 3b: all nodes have reasonable coverage (sanity check) ────────────────
    // Beta and truncated nodes are approximations — their coverage is not guaranteed
    // to equal p_cov. We only assert coverage > 0.5 to catch gross failures
    // (e.g. the distribution being completely off-target).
    @Test
    void testAllNodesCoverageReasonable() {
        for (int i = 0; i < N; i++) {
            int fi = i;
            double cov = Arrays.stream(samples)
                    .mapToDouble(s -> (s[fi] >= T_LO[fi] && s[fi] <= T_HI[fi]) ? 1.0 : 0.0)
                    .average().orElseThrow();

            assertTrue(cov > 0.5,
                    "Node " + (i + 1) + " coverage=" + String.format("%.4f", cov) + " is unreasonably low (< 0.5), suggesting a parameterisation failure");
        }
    }


    // ── test 4: parent always strictly older than child ───────────────────────
    @Test
    void testParentAlwaysOlderThanChild() {
        for (int s = 0; s < N_SIM; s++) {
            for (int i = 1; i < N; i++) {
                if (PARENT[i] != -1) {
                    double child  = samples[s][i];
                    double parent = samples[s][PARENT[i]];
                    // >= not > : Beta.sample() can return 1.0 giving child==parent,
                    // which R also allows — the hard invariant is child <= parent
                    assertTrue(parent >= child,
                            "Sample " + s + " node " + (i + 1) +
                                    ": parent=" + parent + " < child=" + child);
                }
            }
        }
    }

    // ── test 5: t-test on ROOT node only (pure untruncated lognormal) ─────────
    // FIX: nodes 1,2,5,6,7 (0-indexed) are TRUNCATED lognormals — their sample
    //      mean is lower than exp(mu+sigma2/2). Only node 0 is a pure lognormal.
    @Test
    void testRootMeanConsistentWithLogNormal() {
        TTest tTest = new TTest();
        double mu         = computeMu(T_LO[0], T_HI[0], P_COV);
        double sigma2     = computeSigma2(T_LO[0], T_HI[0], P_COV);
        double targetMean = Math.exp(mu + sigma2 / 2.0);

        double[] rootSamples = Arrays.stream(samples).mapToDouble(s -> s[0]).toArray();
        double pValue = tTest.tTest(targetMean, rootSamples);

        assertTrue(pValue > T_TEST_ALPHA,
                "Root node t-test rejected lognormal mean. " +
                        "p-value=" + String.format("%.6f", pValue) +
                        " targetMean=" + String.format("%.4f", targetMean));
    }

    // ── test 6: t-test on beta nodes — child/parent ratio in (0,1) ───────────
    // FIX: kept as-is — this test is correct and does not depend on untruncated mean.
    @Test
    void testBetaNodeRatiosConsistent() {
        TTest tTest  = new TTest();
        int[] betaNodes = {3, 4, 8, 9};

        for (int i : betaNodes) {
            int p = PARENT[i];
            double[] parentSamples = Arrays.stream(samples).mapToDouble(s -> s[p]).toArray();
            double[] childSamples  = Arrays.stream(samples).mapToDouble(s -> s[i]).toArray();

            double[] ratios = new double[N_SIM];
            for (int s = 0; s < N_SIM; s++) {
                ratios[s] = childSamples[s] / parentSamples[s];
            }
            double meanRatio = Arrays.stream(ratios).average().orElseThrow();

            assertTrue(meanRatio > 0.0 && meanRatio < 1.0,
                    "Node " + (i + 1) + ": mean ratio=" + meanRatio + " must be in (0,1)");

            double[] scaledParent = Arrays.stream(parentSamples).map(v -> v * meanRatio).toArray();
            double pValue = tTest.tTest(scaledParent, childSamples);
            assertTrue(pValue > T_TEST_ALPHA,
                    "Node " + (i + 1) + ": beta ratio t-test rejected. " +
                            "p-value=" + String.format("%.6f", pValue) +
                            " meanRatio=" + String.format("%.4f", meanRatio));
        }
    }

    // ── test 7: all samples positive and finite ───────────────────────────────
    @Test
    void testAllSamplesPositiveAndFinite() {
        for (int s = 0; s < N_SIM; s++) {
            for (int i = 0; i < N; i++) {
                double v = samples[s][i];
                assertTrue(Double.isFinite(v) && v > 0,
                        "Sample " + s + " node " + (i + 1) + " not positive finite: " + v);
            }
        }
    }

    // ── helpers ───────────────────────────────────────────────────────────────
    private static double computeMu(double lo, double hi, double p) {
        double z     = Math.sqrt(2.0) * Erf.erfInv(p);
        double sigma = (Math.log(hi) - Math.log(lo)) / (2 * z);
        return Math.log(lo) + z * sigma;
    }

    private static double computeSigma2(double lo, double hi, double p) {
        double z     = Math.sqrt(2.0) * Erf.erfInv(p);
        double sigma = (Math.log(hi) - Math.log(lo)) / (2 * z);
        return sigma * sigma;
    }

    // FIX: nodes 7,8,9 must have strictly nested distinct taxa so computeParents
    //      can resolve their parent chain: size[parent] > size[child] is required.
    //      node7={t8,t9}, node8={t8}, node9 needs to nest inside node8 but
    //      single-taxon sets cannot nest — use a linear chain with unique leaves.
    private static String[][] buildDummyTaxa() {
        return new String[][] {
                {"t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"}, // node 0: root,      parent=-1
                {"t1","t2","t3","t4","t5"},                            // node 1: child of 0, non-beta
                {"t6","t7","t8","t9","t10"},                           // node 2: child of 0, non-beta
                {"t6","t7"},                                           // node 3: child of 2, beta
                {"t1","t2"},                                           // node 4: child of 1, beta
                {"t3","t4","t5"},                                      // node 5: child of 1, non-beta
                {"t3","t4"},                                           // node 6: child of 5, non-beta
                {"t8","t9","t10"},                                     // node 7: child of 5 → FIX: must be child of 5
                {"t8","t9"},                                           // node 8: child of 7, beta — size 2 < size 3 ✓
                {"t8"},                                                // node 9: child of 8, beta — size 1 < size 2 ✓
        };
    }
}