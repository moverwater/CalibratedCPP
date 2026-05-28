import calibratedcpp.lphy.tree.CalibratedAgeDependentCPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.*;

/**
 * Exploratory printer — not an assertion test.
 *
 * Samples trees from CalibratedAgeDependentCPPTree with a fixed rootAge but NO
 * hard calibrations, drawing (lambda, lifetime) from user-specified priors each
 * iteration.  Reports plausible calibration age ranges in three groups:
 *
 *   Group 1 — single 4-tip calibration
 *   Group 2 — disjoint 4-tip + 5-tip calibrations
 *   Group 3 — 3-tip nested inside 6-tip, plus a separate 5-tip calibration
 *
 * Edit the CONFIG block to match your own priors before running.
 */
public class CalibrationRangeExplorerTest {

    // =========================================================================
    // CONFIG — edit these
    // =========================================================================

    private static final int    N_TREES     = 2000;
    private static final int    N_TAXA      = 100;
    private static final double RHO         = 1.0;
    private static final double ROOT_AGE    = 3.0;

    // lambda ~ LogNormal(meanLog, sdLog)
    // median lambda = exp(0.0) = 1.0; splits spread across [0, rootAge] with rho=1 and mu≈0
    private static final double LAMBDA_MEAN_LOG = 0.0;
    private static final double LAMBDA_SD_LOG   = 0.5;

    // lifetime ~ Gamma(shape, scale)  — matched to LogNormal(8, 0.5) for mu
    private static final double LIFETIME_SHAPE  = 6.0;
    private static final double LIFETIME_SCALE  = 9.0;

    private static final long SEED = 42L;

    // =========================================================================
    // Distribution samplers
    // =========================================================================

    static double sampleLogNormal(Random rng, double meanLog, double sdLog) {
        return Math.exp(meanLog + sdLog * rng.nextGaussian());
    }

    /**
     * Gamma sampler via Marsaglia-Tsang (2000).
     * For shape < 1 uses the boost trick: Gamma(shape+1) * U^(1/shape).
     */
    static double sampleGamma(Random rng, double shape, double scale) {
        if (shape < 1.0) {
            return sampleGamma(rng, shape + 1.0, scale)
                    * Math.pow(rng.nextDouble(), 1.0 / shape);
        }
        double d = shape - 1.0 / 3.0;
        double c = 1.0 / Math.sqrt(9.0 * d);
        while (true) {
            double x, v;
            do {
                x = rng.nextGaussian();
                v = 1.0 + c * x;
            } while (v <= 0.0);
            v = v * v * v;
            double u = rng.nextDouble();
            if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v * scale;
            if (Math.log(u) < 0.5 * x * x + d * (1.0 - v + Math.log(v))) return d * v * scale;
        }
    }

    // =========================================================================
    // Tree helpers
    // =========================================================================

    static int leafCount(TimeTreeNode node) {
        if (node.isLeaf()) return 1;
        int c = 0;
        for (TimeTreeNode ch : node.getChildren()) c += leafCount(ch);
        return c;
    }

    /** All internal nodes in the tree that span exactly {@code k} leaves. */
    static List<TimeTreeNode> nodesOfSize(TimeTree tree, int k) {
        List<TimeTreeNode> result = new ArrayList<>();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode nd = tree.getNodeByIndex(i);
            if (!nd.isLeaf() && leafCount(nd) == k) result.add(nd);
        }
        return result;
    }

    /**
     * Returns true if {@code ancestor} is a strict ancestor of {@code descendant}
     * (i.e. {@code ancestor} is reachable by following parent links from {@code descendant},
     * but they are not the same node).
     */
    static boolean isProperAncestor(TimeTreeNode ancestor, TimeTreeNode descendant) {
        TimeTreeNode cur = descendant.getParent();
        while (cur != null) {
            if (cur == ancestor) return true;
            cur = cur.getParent();
        }
        return false;
    }

    /** True if neither node is an ancestor of the other (their leaf sets are disjoint). */
    static boolean areDisjoint(TimeTreeNode a, TimeTreeNode b) {
        return !isProperAncestor(a, b) && !isProperAncestor(b, a);
    }

    /**
     * Finds one pair of internal nodes [nodeK1, nodeK2] that are disjoint,
     * where nodeK1 spans k1 leaves and nodeK2 spans k2 leaves.
     * Returns null if no such pair exists in this tree.
     */
    static TimeTreeNode[] findDisjointPair(TimeTree tree, int k1, int k2) {
        List<TimeTreeNode> listK1 = nodesOfSize(tree, k1);
        List<TimeTreeNode> listK2 = nodesOfSize(tree, k2);
        for (TimeTreeNode a : listK1) {
            for (TimeTreeNode b : listK2) {
                if (areDisjoint(a, b)) return new TimeTreeNode[]{a, b};
            }
        }
        return null;
    }

    /**
     * Finds a triple [outerNode, innerNode, sepNode] where:
     *   - outerNode spans outerK leaves
     *   - innerNode spans innerK leaves and is strictly inside outerNode
     *   - sepNode spans sepK leaves and is disjoint from outerNode
     * Returns null if no such triple exists in this tree.
     */
    static TimeTreeNode[] findNestedPlusDisjoint(TimeTree tree, int outerK, int innerK, int sepK) {
        List<TimeTreeNode> outerList = nodesOfSize(tree, outerK);
        List<TimeTreeNode> innerList = nodesOfSize(tree, innerK);
        List<TimeTreeNode> sepList   = nodesOfSize(tree, sepK);

        for (TimeTreeNode outer : outerList) {
            for (TimeTreeNode inner : innerList) {
                if (!isProperAncestor(outer, inner)) continue;
                for (TimeTreeNode sep : sepList) {
                    if (areDisjoint(outer, sep)) {
                        return new TimeTreeNode[]{outer, inner, sep};
                    }
                }
            }
        }
        return null;
    }

    /** Linear-interpolation percentile on a sorted array. */
    static double percentile(double[] sorted, double p) {
        if (sorted.length == 0) return Double.NaN;
        double idx = p / 100.0 * (sorted.length - 1);
        int lo = (int) idx;
        int hi = Math.min(lo + 1, sorted.length - 1);
        return sorted[lo] + (idx - lo) * (sorted[hi] - sorted[lo]);
    }

    static void printRow(String label, List<Double> ages) {
        if (ages.isEmpty()) {
            System.out.printf("  %-30s  %6d  (no data)%n", label, 0);
            return;
        }
        double[] sorted = ages.stream().mapToDouble(Double::doubleValue).sorted().toArray();
        System.out.printf("  %-30s  %6d  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f%n",
                label, sorted.length,
                percentile(sorted, 5),
                percentile(sorted, 25),
                percentile(sorted, 50),
                percentile(sorted, 75),
                percentile(sorted, 95));
    }

    // =========================================================================
    // Test
    // =========================================================================

    @Test
    void printCalibrationRanges() {
        Random rng = new Random(SEED);
        int skipped = 0;

        // Group 1: single 4-tip calibration
        List<Double> g1_4tip = new ArrayList<>();

        // Group 2: disjoint 4-tip + 5-tip
        List<Double> g2_4tip = new ArrayList<>();
        List<Double> g2_5tip = new ArrayList<>();

        // Group 3: 3-tip nested inside 6-tip, plus separate 5-tip
        List<Double> g3_6tip = new ArrayList<>();
        List<Double> g3_3tip = new ArrayList<>();
        List<Double> g3_5tip = new ArrayList<>();

        for (int i = 0; i < N_TREES; i++) {
            double lambda   = sampleLogNormal(rng, LAMBDA_MEAN_LOG, LAMBDA_SD_LOG);
            double lifetime = sampleGamma(rng, LIFETIME_SHAPE, LIFETIME_SCALE);

            CalibratedAgeDependentCPPTree gen = new CalibratedAgeDependentCPPTree(
                    new Value<>("", lambda),
                    new Value<>("", RHO),
                    new Value<>("", N_TAXA),
                    new Value<>("", lifetime),
                    null, null, null,
                    new Value<>("", ROOT_AGE));

            TimeTree tree;
            try {
                tree = gen.sample().value();
            } catch (Exception e) {
                skipped++;
                continue;
            }

            // Group 1: any 4-tip internal node
            List<TimeTreeNode> nodes4 = nodesOfSize(tree, 4);
            if (!nodes4.isEmpty()) {
                g1_4tip.add(nodes4.get(0).getAge());
            }

            // Group 2: a disjoint 4-tip + 5-tip pair
            TimeTreeNode[] pair = findDisjointPair(tree, 4, 5);
            if (pair != null) {
                g2_4tip.add(pair[0].getAge());
                g2_5tip.add(pair[1].getAge());
            }

            // Group 3: 3-tip nested in 6-tip, plus disjoint 5-tip
            TimeTreeNode[] triple = findNestedPlusDisjoint(tree, 6, 3, 5);
            if (triple != null) {
                g3_6tip.add(triple[0].getAge());
                g3_3tip.add(triple[1].getAge());
                g3_5tip.add(triple[2].getAge());
            }
        }

        // ── Header ────────────────────────────────────────────────────────────
        System.out.println();
        System.out.println("=== Calibration Range Explorer ===");
        System.out.printf("  lambda   ~ LogNormal(meanLog=%.2f, sdLog=%.2f)   [median=%.3f]%n",
                LAMBDA_MEAN_LOG, LAMBDA_SD_LOG, Math.exp(LAMBDA_MEAN_LOG));
        System.out.printf("  lifetime ~ Gamma(shape=%.2f, scale=%.2f)   [mean=%.1f]%n",
                LIFETIME_SHAPE, LIFETIME_SCALE, LIFETIME_SHAPE * LIFETIME_SCALE);
        System.out.printf("  rootAge=%.1f   n=%d   rho=%.2f   N_TREES=%d   skipped=%d%n%n",
                ROOT_AGE, N_TAXA, RHO, N_TREES, skipped);

        String header = String.format("  %-30s  %6s  %7s  %7s  %7s  %7s  %7s",
                "Calibration", "n", "5th", "25th", "Median", "75th", "95th");
        String rule = "  " + "-".repeat(header.length() - 2);

        // ── Group 1 ───────────────────────────────────────────────────────────
        System.out.println("--- Group 1: Single 4-tip calibration ---");
        System.out.println(header);
        System.out.println(rule);
        printRow("4-tip MRCA age", g1_4tip);
        System.out.println();

        // ── Group 2 ───────────────────────────────────────────────────────────
        System.out.println("--- Group 2: Disjoint 4-tip + 5-tip calibrations ---");
        System.out.println(header);
        System.out.println(rule);
        printRow("4-tip MRCA age", g2_4tip);
        printRow("5-tip MRCA age", g2_5tip);
        System.out.println();

        // ── Group 3 ───────────────────────────────────────────────────────────
        System.out.println("--- Group 3: 3-tip nested in 6-tip + separate 5-tip ---");
        System.out.println(header);
        System.out.println(rule);
        printRow("6-tip outer MRCA age", g3_6tip);
        printRow("3-tip inner MRCA age (< 6-tip)", g3_3tip);
        printRow("5-tip separate MRCA age", g3_5tip);
        System.out.println();
    }
}
