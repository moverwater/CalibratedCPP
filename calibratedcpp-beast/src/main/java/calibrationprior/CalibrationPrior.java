package calibrationprior;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import calibration.CalibrationForest;
import calibration.CalibrationNode;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;

import java.util.*;

/**
 * @author Marcus Overwater
 */

@Description("A topology consistent prior on clade calibration ages.")
public class CalibrationPrior extends Distribution {

    public Input<TreeInterface> treeInput =
            new Input<>("tree", "Tree to calibrate", Input.Validate.REQUIRED);

    public Input<List<CalibrationCladePrior>> cladesInput =
            new Input<>("calibration", "List of calibration clades", new ArrayList<>());

    private List<CalibrationNode> calibrationNodes = new ArrayList<>();

    @Override
    public void initAndValidate() {
        TreeInterface tree = treeInput.get();
        List<CalibrationCladePrior> cladePriors = cladesInput.get();
        if (tree == null) throw new IllegalArgumentException("Tree is null");
        if (cladePriors.isEmpty()) {
            throw new IllegalArgumentException("No clades provided");
        }

        // === Step 2: build inclusion hierarchy ===
        CalibrationForest calibrationForest = new CalibrationForest(cladePriors);
        calibrationNodes = calibrationForest.getAllNodes();

        for (CalibrationNode node : calibrationNodes) {
            if (node.parent != null) {
                node.getCalibrationCladePrior().isOverlapEdge = node.parent.getCalibrationCladePrior().getLower() < node.getCalibrationCladePrior().getUpper();
            }
        }

        // === Step 3: compute log-moments ===
        for (CalibrationNode n : calibrationNodes) computeLogTargets(n);

        // === Step 4: partition into overlap components and fit ===
        List<CalibrationComponent> comps = CalibrationComponent.partition(calibrationForest);
        for (CalibrationComponent comp : comps) fitComponent(comp);
        // done
    }

    private void fitComponent(CalibrationComponent comp) {
        CalibrationNode root = comp.getRoot();
        List<CalibrationNode> members = comp.getMembers();

        // Map each node to a column index in the least-squares system
        Map<CalibrationNode, Integer> colIndex = new LinkedHashMap<>();
        int idx = 0;
        for (CalibrationNode n : members) {
            if (n != root) colIndex.put(n, idx++);
        }

        int K = colIndex.size();
        if (K == 0) return; // only root in component, nothing to fit

        RealMatrix A = new Array2DRowRealMatrix(K, K);
        double[] bMu = new double[K];
        double[] bVar = new double[K];

        // Build linear system: log(child) - log(root) = sum of edges along path
        for (CalibrationNode child : colIndex.keySet()) {
            int row = colIndex.get(child);

            List<CalibrationNode> path = pathFromRoot(root, child);
            for (CalibrationNode step : path) {
                if (step == root) continue;
                Integer col = colIndex.get(step);
                if (col != null) {
                    A.addToEntry(row, col, 1.0);
                }
            }

            bMu[row] = child.getCalibrationCladePrior().mu - root.getCalibrationCladePrior().mu;
            bVar[row] = child.getCalibrationCladePrior().sigma2 - root.getCalibrationCladePrior().sigma2;
        }

        // Solve least-squares for mean and variance edges
        RealVector mHat = solveQR(A, new ArrayRealVector(bMu));
        RealVector vHat = solveQR(A, new ArrayRealVector(bVar));

        for (int i = 0; i < K; i++) {
            if (vHat.getEntry(i) <= 0) vHat.setEntry(i, 1e-8);
        }

        // Assign edges and Beta parameters
        for (Map.Entry<CalibrationNode, Integer> e : colIndex.entrySet()) {
            CalibrationNode child = e.getKey();
            int j = e.getValue();

            child.getCalibrationCladePrior().mEdge = mHat.getEntry(j);
            child.getCalibrationCladePrior().vEdge = vHat.getEntry(j);

            double[] ab = invertLogMomentsToBetaParams(child.getCalibrationCladePrior().mEdge, child.getCalibrationCladePrior().vEdge);
            child.getCalibrationCladePrior().alpha = ab[0];
            child.getCalibrationCladePrior().beta = ab[1];
        }
    }

    // ------------------------------------------------------------------
    // 3. Fit each component

    private RealVector solveQR(RealMatrix A, RealVector b) {
        try {
            return new QRDecomposition(A).getSolver().solve(b);
        } catch (Exception e) {
            return new SingularValueDecomposition(A).getSolver().solve(b);
        }
    }

    private List<CalibrationNode> pathFromRoot(CalibrationNode root, CalibrationNode node) {
        LinkedList<CalibrationNode> path = new LinkedList<>();
        CalibrationNode cur = node;
        while (cur != null && cur != root) {
            path.addFirst(cur);
            cur = cur.parent;
        }
        path.addFirst(root);
        return path;
    }

    // ------------------------------------------------------------------
    // log-moment target from calibration bounds
    private void computeLogTargets(CalibrationNode n) {
        double tLo = n.getCalibrationCladePrior().getLower();
        double tHi = n.getCalibrationCladePrior().getUpper();
        double p = n.getCalibrationCladePrior().getCoverage();

        // Replaced NormalDistribution.inverseCumulativeProbability with Erf.erfInv
        // This calculates the z-score for the given coverage probability
        double z = Math.sqrt(2.0) * Erf.erfInv(p);

        double sigma = (Math.log(tHi) - Math.log(tLo)) / (2 * z);
        n.getCalibrationCladePrior().sigma2 = sigma * sigma;
        n.getCalibrationCladePrior().mu = Math.log(tLo) + z * sigma;
    }

    // ------------------------------------------------------------------
    // Utility math
    private double[] invertLogMomentsToBetaParams(double m, double v) {
        double Ey = Math.exp(m + v / 2);
        double Vy = Math.exp(2 * m + v) * (Math.exp(v) - 1);
        if (Vy <= 0) Vy = 1e-8;
        Ey = Math.min(1 - 1e-8, Math.max(1e-8, Ey));

        double common = Ey * (1 - Ey) / Vy - 1;
        double a = Math.max(1e-3, Ey * common);
        double b = Math.max(1e-3, (1 - Ey) * common);

        for (int iter = 0; iter < 40; iter++) {
            double psiA = Gamma.digamma(a);
            double psiAB = Gamma.digamma(a + b);
            double f1 = psiA - psiAB - m;
            double f2 = Gamma.trigamma(a) - Gamma.trigamma(a + b) - v;
            if (Math.abs(f1) < 1e-8 && Math.abs(f2) < 1e-8) break;
            a -= 0.5 * f1;
            b -= 0.5 * f2;
            a = Math.max(a, 1e-6);
            b = Math.max(b, 1e-6);
        }
        return new double[]{a, b};
    }

    private void collectLeafTaxa(Node node, Set<String> leafIDs) {
        if (node.isLeaf()) {
            leafIDs.add(node.getID());
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                collectLeafTaxa(node.getChild(i), leafIDs);
            }
        }
    }

    // ------------------------------------------------------------------
    // BEAST interface
    @Override
    public double calculateLogP() {
        TreeInterface tree = treeInput.get();
        double logP = 0;
        for (CalibrationNode n : calibrationNodes) {

            // Check for monophyly
            Set<String> leafIDs = new HashSet<>();
            Node treeNode = n.getCommonAncestor(tree);
            collectLeafTaxa(treeNode, leafIDs);
            if (!leafIDs.equals(n.getCalibrationClade().getTaxa().getTaxaNames())) {
                return Double.NEGATIVE_INFINITY; // clade is not monophyletic!
            }

            double t = treeNode.getHeight();
            if (n.parent == null) {
                // root lognormal
                double lp = logNormalLogPdf(t, n.getCalibrationCladePrior().mu, Math.sqrt(n.getCalibrationCladePrior().sigma2));
                logP += lp;
            } else {
                Node pNode = n.parent.getCommonAncestor(tree);
                double tp = pNode.getHeight();
                if (n.getCalibrationCladePrior().isOverlapEdge) {
                    // overlapping → Beta
                    double r = t / tp;
                    if (r <= 0 || r >= 1) return Double.NEGATIVE_INFINITY;
                    logP += betaLogPdf(r, n.getCalibrationCladePrior().alpha, n.getCalibrationCladePrior().beta) - Math.log(tp); // Jacobian for conversion of the joint density of tp and r to tp and t
                } else {
                    // non-overlapping → truncated lognormal
                    double lp = logNormalLogPdf(t, n.getCalibrationCladePrior().mu, Math.sqrt(n.getCalibrationCladePrior().sigma2));
                    double lcdf = logNormalLogCdf(tp, n.getCalibrationCladePrior().mu, Math.sqrt(n.getCalibrationCladePrior().sigma2));
                    if (Double.isInfinite(lcdf)) return Double.NEGATIVE_INFINITY;
                    logP += lp - lcdf;
                }
            }
        }
        return logP;
    }

    private double logNormalLogPdf(double x, double mu, double sigma) {
        if (x <= 0) return Double.NEGATIVE_INFINITY;
        double z = (Math.log(x) - mu) / sigma;
        return -Math.log(x * sigma * Math.sqrt(2 * Math.PI)) - 0.5 * z * z;
    }

    private double logNormalLogCdf(double x, double mu, double sigma) {
        if (x <= 0) return Double.NEGATIVE_INFINITY;

        // 1. Convert to Z-score
        double z = (Math.log(x) - mu) / sigma;

        // 2. Calculate Cumulative Probability using Erf
        // Formula: 0.5 * (1 + erf(z / sqrt(2)))
        double p = 0.5 * (1.0 + Erf.erf(z / Math.sqrt(2.0)));

        return Math.log(p);
    }

    private double betaLogPdf(double x, double a, double b) {
        if (x <= 0 || x >= 1) return Double.NEGATIVE_INFINITY;
        return (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x)
                - (Gamma.logGamma(a) + Gamma.logGamma(b) - Gamma.logGamma(a + b));
    }

    public List<CalibrationCladePrior> getCalibrationCladePriors() {
        return cladesInput.get();
    }

    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    protected boolean requiresRecalculation() {
        // Do any updates
        for (CalibrationNode calibration : calibrationNodes) {
            Node mrca = calibration.getCommonAncestor(treeInput.get());
            if (calibration.getCalibrationCladePrior().getAge() != null) {
                calibration.getCalibrationClade().getAge().setValue(mrca.getHeight());
            }
        }
        return true;
    }

    // ------------------------------------------------------------------
    // 2. Decompose forest into subtrees of clades that have intervals overlapping with their parents
    public static class CalibrationComponent {
        private final List<CalibrationNode> members = new ArrayList<>();
        private CalibrationNode root; // exactly one root per component

        public CalibrationComponent() {
        }

        public static List<CalibrationComponent> partition(CalibrationForest calibrationForest) {
            List<CalibrationNode> all = calibrationForest.getAllNodes();
            List<CalibrationComponent> components = new ArrayList<>();
            Set<CalibrationNode> visited = new HashSet<>();

            for (CalibrationNode start : all) {
                if (visited.contains(start)) continue;

                CalibrationComponent comp = new CalibrationComponent();
                Deque<CalibrationNode> stack = new ArrayDeque<>();
                stack.push(start);
                visited.add(start);

                // --- DFS over overlapping edges ---
                while (!stack.isEmpty()) {
                    CalibrationNode cur = stack.pop();
                    comp.add(cur);

                    for (CalibrationNode other : all) {
                        if (cur == other || visited.contains(other)) continue;
                        if (overlaps(cur, other)) {
                            stack.push(other);
                            visited.add(other);
                        }
                    }
                }

                // --- Identify the highest node in the component ---
                comp.root = comp.findComponentRoot();
                components.add(comp);
            }

            return components;
        }

        private static boolean overlaps(CalibrationNode a, CalibrationNode b) {
            if (a.parent == b)
                return b.getCalibrationCladePrior().getLower() < a.getCalibrationCladePrior().getUpper(); // child inside parent's range
            if (b.parent == a)
                return a.getCalibrationCladePrior().getLower() < b.getCalibrationCladePrior().getUpper();
            return false;
        }

        public void add(CalibrationNode n) {
            members.add(n);
        }

        public CalibrationNode getRoot() {
            return root;
        }

        public List<CalibrationNode> getMembers() {
            return members;
        }

        private CalibrationNode findComponentRoot() {
            CalibrationNode rootCandidate = null;
            for (CalibrationNode n : members) {
                // A root is one whose parent is either null or not in the component,
                // or one whose parent does not overlap with it.
                if (n.parent == null || !members.contains(n.parent) || !overlaps(n, n.parent)) {
                    if (rootCandidate != null)
                        throw new IllegalStateException("Multiple roots found in component!");
                    rootCandidate = n;
                }
            }
            if (rootCandidate == null)
                throw new IllegalStateException("No root found for component!");
            return rootCandidate;
        }
    }
}
