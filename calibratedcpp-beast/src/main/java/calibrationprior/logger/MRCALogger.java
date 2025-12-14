package calibrationprior.logger;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;

import java.io.PrintStream;
import java.util.*;

/**
 * @author Marcus Overwater
 */

@Description("Logs the MRCA ages of calibration clades")
public class MRCALogger extends CalculationNode implements Loggable {

    public Input<TreeInterface> treeInput =
            new Input<>("tree", "The tree whose MRCA ages should be logged", Input.Validate.REQUIRED);

    public Input<CalibrationPrior> cladesInput =
            new Input<>("calibrations", "The calibration clades whose MRCA ages should be logged", Input.Validate.REQUIRED);

    private TreeInterface tree;
    private List<CalibrationCladePrior> clades;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        clades = cladesInput.get().getCalibrationCladePriors();
        if (tree == null)
            throw new IllegalArgumentException("Tree is null");
        if (clades.isEmpty())
            throw new IllegalArgumentException("No clades provided to MRCALogger");
    }

    @Override
    public void init(PrintStream out) {
        for (CalibrationCladePrior c : clades) {
            out.print(getID() + "." + c.getTaxa().getID() + "\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (CalibrationCladePrior c : clades) {
            Node mrca = findMRCA(tree, c.getTaxa());
            out.print(mrca.getHeight() + "\t");
        }
    }

    @Override
    public void close(PrintStream out) { }

    // ----------------------- Helper methods -----------------------

    private Node findMRCA(TreeInterface tree, TaxonSet taxonSet) {
        // 1. Optimization: Index all leaf nodes in the tree by their ID.
        // This allows O(1) lookup instead of iterating through the whole tree for every taxon.
        Map<String, Node> leafMap = new HashMap<>();
        for (Node n : tree.getExternalNodes()) {
            leafMap.put(n.getID(), n);
        }

        List<Node> leafNodes = new ArrayList<>();
        for (Taxon t : taxonSet.getTaxonSet()) {
            // O(1) lookup
            Node n = leafMap.get(t.getID());

            // Safety check: ensure the taxon actually exists in the tree
            if (n != null) {
                leafNodes.add(n);
            }
        }

        if (leafNodes.isEmpty()) return null;

        // 2. Reduce logic remains the same (Iterative MRCA)
        Node mrca = leafNodes.get(0);
        for (int i = 1; i < leafNodes.size(); i++) {
            mrca = getMRCA(mrca, leafNodes.get(i));
            // Optimization: If we hit the root, we can stop immediately.
            assert mrca != null;
            if (mrca.isRoot()) return mrca;
        }
        return mrca;
    }

    private Node getMRCA(Node a, Node b) {
        Set<Node> ancestors = new HashSet<>();
        Node cur = a;
        while (cur != null) {
            ancestors.add(cur);
            cur = cur.getParent();
        }
        cur = b;
        while (cur != null) {
            if (ancestors.contains(cur)) return cur;
            cur = cur.getParent();
        }
        return null;
    }
}
