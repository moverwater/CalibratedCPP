package calibration;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import calibrationprior.CalibrationCladePrior;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Marcus Overwater
 */

@Description("A calibration node is a node in a calibration forest.")
public class CalibrationNode extends BEASTObject {
    public CalibrationClade calibration;
    public TaxonSet taxa;

    public CalibrationNode parent;
    public List<CalibrationNode> children = new ArrayList<>();

    public boolean isRoot;        // true if start of independent subtree

    public CalibrationNode(CalibrationClade calibration) {
        this.calibration = calibration;
        this.taxa = calibration.getTaxa();
    }

    public static CalibrationNode getByTaxa(List<CalibrationNode> nodes, TaxonSet taxa) {
        Set<String> targetTaxa = new HashSet<>();
        for (Taxon t : taxa.getTaxonSet()) targetTaxa.add(t.getID());

        for (CalibrationNode node : nodes) {
            Set<String> nodeTaxa = new HashSet<>();
            for (Taxon t : node.taxa.getTaxonSet()) nodeTaxa.add(t.getID());
            if (nodeTaxa.equals(targetTaxa)) return node;
        }
        return null;
    }

    public CalibrationClade getCalibrationClade() {
        return calibration;
    }

    public CalibrationCladePrior getCalibrationCladePrior() {
        if (!(calibration instanceof CalibrationCladePrior)) {
            throw new IllegalStateException("Calibration node does not hold a CalibrationCladePrior");
        }
        return (CalibrationCladePrior) calibration;
    }

    // --- Utility ---
    public Node getCommonAncestor(TreeInterface tree) {
        List<Node> nodes = new ArrayList<>();
        for (Taxon t : this.taxa.getTaxonSet()) {
            int idx = tree.getTaxonset().getTaxonIndex(t.getID());
            if (idx < 0) throw new IllegalArgumentException("Taxon " + t.getID() + " not in tree");
            nodes.add(tree.getNode(idx));
        }
        if (nodes.isEmpty()) return null;

        Node a = nodes.get(0);
        for (int i = 1; i < nodes.size(); i++) {
            a = getCommonAncestor(a, nodes.get(i));
        }
        return a;
    }

    private Node getCommonAncestor(Node a, Node b) {
        if (a == b && a.getParent() == a) {
            throw new IllegalStateException("Node " + a + " is its own parent!");
        }
        Set<Node> ancA = new HashSet<>();
        while (a != null) { ancA.add(a); a = a.getParent(); }
        while (b != null) {
            if (ancA.contains(b)) return b;
            b = b.getParent();
        }
        return null;
    }

    @Override
    public String toString() {
        String id = taxa != null ? taxa.getID() : "noID";
        return String.format("CalibrationNode[%s]", id);
    }

    @Override
    public void initAndValidate() {

    }
}
