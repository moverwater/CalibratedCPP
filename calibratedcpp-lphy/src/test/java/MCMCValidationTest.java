import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.tree.CalibratedAgeDependentCPPTree;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.Tree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.WrappedJEBLTimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Java equivalent of the R validation scripts:
 *   lphy_added_fixed_stem_validation_100.R
 *   lphy_added_fixed_root_validation_100.R
 *
 * Each test runs two stages:
 *   Stage 1 — MCMC trees vs CalibratedCPPTree forward simulation (Welch's t-test).
 *   Stage 2 — CalibratedCPPTree vs CalibratedAgeDependentCPPTree (KS test).
 *
 * Tests are skipped silently when the MCMC file is absent.
 */
public class MCMCValidationTest {

    private static final int N_SIM = 1000;
    private static final int THIN  = 5;    // keep every 5th MCMC tree to reduce autocorrelation

    // =========================================================================
    // NEXUS tree loading
    // =========================================================================

    /**
     * Reads all trees from a BEAST NEXUS file and keeps every {@code thin}-th one.
     * Thinning reduces autocorrelation so the nominal SE in Welch's t-test is a
     * better approximation of the true sampling uncertainty.
     */
    static TimeTree[] readNexusTrees(String path, int thin) throws IOException, ImportException {
        try (BufferedReader reader = Files.newBufferedReader(Path.of(path), StandardCharsets.UTF_8)) {
            NexusImporter importer = new NexusImporter(reader, false, 0);
            List<Tree> raw = importer.importTrees();
            int kept = (raw.size() + thin - 1) / thin;
            TimeTree[] trees = new TimeTree[kept];
            int idx = 0;
            for (int i = 0; i < raw.size(); i += thin)
                trees[idx++] = WrappedJEBLTimeTreeNode.Utils.convert(raw.get(i));
            return trees;
        }
    }

    // =========================================================================
    // Stat collection for pre-generated tree arrays
    // =========================================================================

    static CalibratedCPPDistributionTest.Stats computeStats(TimeTree[] trees) {
        int n = trees.length;
        double[] len = new double[n], root = new double[n],
                 mean = new double[n], maxNR = new double[n],
                 gam  = new double[n], col  = new double[n],
                 cher = new double[n], bet  = new double[n];

        for (int i = 0; i < n; i++) {
            TimeTree t = trees[i];
            len[i]  = CalibratedCPPDistributionTest.treeLength(t);
            gam[i]  = CalibratedCPPDistributionTest.gammaStatistic(t);
            col[i]  = CalibratedCPPDistributionTest.collessIndex(t);
            cher[i] = CalibratedCPPDistributionTest.numCherries(t);
            double b = CalibratedCPPDistributionTest.betaStatistic(t);
            bet[i]  = Double.isFinite(b) ? b : Double.NaN;
            double[] ages = CalibratedCPPDistributionTest.sortedInternalAges(t);
            root[i]  = ages[ages.length - 1];
            int nr   = ages.length - 1;
            mean[i]  = nr > 0 ? Arrays.stream(ages, 0, nr).average().orElse(0) : 0;
            maxNR[i] = nr > 0 ? ages[nr - 1] : 0;
        }
        return new CalibratedCPPDistributionTest.Stats(len, root, mean, maxNR, gam, col, cher, bet);
    }

    // =========================================================================
    // Welch's t-test
    // =========================================================================

    static double welchT(double[] a, double[] b) {
        double[] fa = Arrays.stream(a).filter(Double::isFinite).toArray();
        double[] fb = Arrays.stream(b).filter(Double::isFinite).toArray();
        int na = fa.length, nb = fb.length;
        double ma = Arrays.stream(fa).average().orElse(0);
        double mb = Arrays.stream(fb).average().orElse(0);
        double va = 0, vb = 0;
        for (double x : fa) va += (x - ma) * (x - ma);
        for (double x : fb) vb += (x - mb) * (x - mb);
        va /= (na - 1);
        vb /= (nb - 1);
        double se = Math.sqrt(va / na + vb / nb);
        // Degenerate: both arrays are constant (e.g. beta hits boundary for large trees).
        if (se < 1e-10) return Math.abs(ma - mb) < 1e-10 ? 0.0 : Double.POSITIVE_INFINITY;
        return (ma - mb) / se;
    }

    /**
     * Asserts |t| < 3.0 (α ≈ 0.003, large-sample normal approximation).
     * The threshold is wider than 1.96 (α = 0.05) to account for MCMC autocorrelation
     * inflating |t| even when distributions match, while still reliably detecting
     * real distributional biases.
     */
    static void assertTTest(double[] a, double[] b, String label) {
        double t = welchT(a, b);
        assertTrue(Math.abs(t) < 3.0,
                String.format("t-test FAILED '%s': t = %.4f (|t| must be < 3.0)", label, t));
    }

    // =========================================================================
    // Shared two-stage validation runner
    // =========================================================================

    /**
     * Runs both validation stages for one scenario.
     *
     * @param label       short name for error messages (e.g. "fixed-stem")
     * @param mcmcPath    absolute path to the BEAST .trees file
     * @param cppSampler  supplier that draws one tree from CalibratedCPPTree
     * @param adSampler   supplier that draws one tree from CalibratedAgeDependentCPPTree
     */
    static void runValidation(String label,
                              String mcmcPath,
                              Supplier<TimeTree> cppSampler,
                              Supplier<TimeTree> adSampler)
            throws IOException, ImportException {

        assumeTrue(new File(mcmcPath).exists(),
                "MCMC file not found (" + mcmcPath + ") — skipping " + label + " validation");

        // ── Load MCMC trees ────────────────────────────────────────────────────
        CalibratedCPPDistributionTest.Stats mcmcStats =
                computeStats(readNexusTrees(mcmcPath, THIN));

        // ── Simulate CalibratedCPPTree trees ──────────────────────────────────
        CalibratedCPPDistributionTest.Stats cppStats =
                CalibratedCPPDistributionTest.collect(N_SIM, cppSampler);

        // ── Stage 1: MCMC vs CalibratedCPPTree (Welch's t-test) ───────────────
        assertTTest(mcmcStats.rootAge(),  cppStats.rootAge(),  label + " MCMC vs CPP: root age");
        assertTTest(mcmcStats.length(),   cppStats.length(),   label + " MCMC vs CPP: tree length");
        assertTTest(mcmcStats.gamma(),    cppStats.gamma(),    label + " MCMC vs CPP: gamma");
        assertTTest(mcmcStats.colless(),  cppStats.colless(),  label + " MCMC vs CPP: Colless");
        assertTTest(mcmcStats.cherries(), cppStats.cherries(), label + " MCMC vs CPP: cherries");
        assertTTest(mcmcStats.beta(),     cppStats.beta(),     label + " MCMC vs CPP: beta");

        // ── Simulate CalibratedAgeDependentCPPTree trees ──────────────────────
        CalibratedCPPDistributionTest.Stats adStats =
                CalibratedCPPDistributionTest.collect(N_SIM, adSampler);

        // ── Stage 2: CalibratedCPPTree vs AgeDependentCPP (KS test) ──────────
        CalibratedCPPDistributionTest.assertKS(cppStats.rootAge(),  adStats.rootAge(),  label + " CPP vs AD: root age");
        CalibratedCPPDistributionTest.assertKS(cppStats.length(),   adStats.length(),   label + " CPP vs AD: tree length");
        CalibratedCPPDistributionTest.assertKS(cppStats.gamma(),    adStats.gamma(),    label + " CPP vs AD: gamma");
        CalibratedCPPDistributionTest.assertKS(cppStats.colless(),  adStats.colless(),  label + " CPP vs AD: Colless");
        CalibratedCPPDistributionTest.assertKS(cppStats.cherries(), adStats.cherries(), label + " CPP vs AD: cherries");
        CalibratedCPPDistributionTest.assertKS(cppStats.beta(),     adStats.beta(),     label + " CPP vs AD: beta");
    }

    // =========================================================================
    // fixed_stem_validation  (λ=2, μ=1, ρ=0.1, stemAge=3, 4 calibrations)
    // =========================================================================

    private static final String STEM_MCMC_PATH ="./mcmcResults/fixStemMCMC.trees";

    private static Calibration[] buildStemCalibrations() {
        String[] clade1 = {"leaf_1", "leaf_2"};
        String[] clade2 = new String[7];
        for (int i = 0; i < 7; i++) clade2[i] = "leaf_" + (45 + i);   // leaf_45..leaf_51
        String[] clade3 = new String[11];
        for (int i = 0; i < 11; i++) clade3[i] = "leaf_" + (45 + i);  // leaf_45..leaf_55
        String[] clade4 = {"leaf_87", "leaf_88", "leaf_89", "leaf_90"};
        return new Calibration[]{
                new Calibration(clade1, 1.0),
                new Calibration(clade2, 1.2),
                new Calibration(clade3, 1.5),
                new Calibration(clade4, 2.0)
        };
    }

    @Test
    void fixedStemValidation() throws IOException, ImportException {
        final double lambda = 2.0, mu = 1.0, rho = 0.1, lifetime = 1.0 / mu;
        final int n = 100;
        final double stemAge = 3.0;
        Calibration[] cals = buildStemCalibrations();

        CalibratedCPPTree cppGen = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(cals)),
                null, new Value<>("", stemAge));

        CalibratedAgeDependentCPPTree adGen = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(cals)),
                null, new Value<>("", stemAge));

        runValidation("fixed-stem", STEM_MCMC_PATH,
                () -> cppGen.sample().value(),
                () -> adGen.sample().value());
    }

    // =========================================================================
    // fixed_root_validation  (λ=0.1, μ=0, ρ=1, root at 10, clade at 3)
    // =========================================================================

    private static final String ROOT_MCMC_PATH ="./mcmcResults/fixRootMCMC.trees";

    private static Calibration[] buildRootCalibrations() {
        // Root calibration: all 100 taxa named "0".."99"
        String[] allTaxa = new String[100];
        for (int i = 0; i < 100; i++) allTaxa[i] = String.valueOf(i);
        // Clade calibration: "45".."55" (11 taxa)
        String[] clade = new String[11];
        for (int i = 0; i < 11; i++) clade[i] = String.valueOf(45 + i);
        return new Calibration[]{
                new Calibration(allTaxa, 10.0),
                new Calibration(clade,    3.0)
        };
    }

    @Test
    void fixedRootValidation() throws IOException, ImportException {
        final double lambda = 0.1, mu = 0.0, rho = 1.0;
        // mu=0 → pure Yule; lifetime=1/mu=∞ → effectiveDeathRate=0 in age-dependent model
        final double lifetime = Double.POSITIVE_INFINITY;
        final int n = 100;
        Calibration[] cals = buildRootCalibrations();

        CalibratedCPPTree cppGen = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(cals)),
                null, null);   // no stem age — root conditioned

        CalibratedAgeDependentCPPTree adGen = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(cals)),
                null, null);   // no stem age — root conditioned

        runValidation("fixed-root", ROOT_MCMC_PATH,
                () -> cppGen.sample().value(),
                () -> adGen.sample().value());
    }
}
