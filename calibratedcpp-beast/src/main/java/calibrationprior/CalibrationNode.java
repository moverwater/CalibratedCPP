package calibrationprior;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Marcus Overwater
 */

@Description("A calibration node is a node in a calibration forest.")
public class CalibrationNode extends BEASTObject {
    public CalibrationClade clade;
    public TaxonSet taxa;

    public CalibrationNode parent;
    public List<CalibrationNode> children = new ArrayList<>();

    public boolean isRoot;        // true if start of independent subtree

    // constructor used in CalibrationPrior
    public CalibrationNode(CalibrationClade clade) {
        this.clade = clade;
        this.taxa = clade.getTaxa();
    }

    // constructor used for CalibratedCPP
    public CalibrationNode(TaxonSet taxa) {
        this.taxa = taxa;
        this.clade = null;
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
        return clade;
    }

    public boolean hasCalibration() {
        return clade != null;
    }

    // --- Utility ---
    public Node getCommonAncestor(TreeInterface tree) {
        List<Node> nodes = new ArrayList<>();
        for (Taxon t : this.taxa.getTaxonSet()) {
            int idx = tree.getTaxonset().getTaxonIndex(t.toString());
            nodes.add(tree.getNode(idx));
        }
        if (nodes.isEmpty()) return null;
        Node a = nodes.get(0);
        for (int i = 1; i < nodes.size(); i++) a = getCommonAncestor(a, nodes.get(i));
        return a;
    }

    private Node getCommonAncestor(Node a, Node b) {
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
        if (hasCalibration()) {
            return String.format("CalibrationNode[%s, [%.3f, %.3f]]", id, clade.getLower(), clade.getUpper());
        } else {
            return String.format("CalibrationNode[%s]", id);
        }
    }

    @Override
    public void initAndValidate() {

    }
}
