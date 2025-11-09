package calibrationprior.logger;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
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
        List<String> taxaNames = tree.getTaxonset().asStringList();
        List<Node> leafNodes = new ArrayList<>();

        for (String taxon : taxonSet.asStringList()) {
            int taxonIndex = taxaNames.indexOf(taxon);
            if (taxonIndex < 0) {
                throw new IllegalArgumentException("Taxon " + taxon + " not found in tree taxa: " + taxaNames);
            }
            Node node = tree.getNode(taxonIndex);
            leafNodes.add(node);
        }

        Node mrca = leafNodes.get(0);
        for (int i = 1; i < leafNodes.size(); i++) {
            mrca = getMRCA(mrca, leafNodes.get(i));
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
