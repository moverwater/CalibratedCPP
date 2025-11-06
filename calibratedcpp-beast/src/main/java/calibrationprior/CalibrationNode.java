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

@Description("A calibration node is a node in a calibration forest for optimising a joint prior on ages of monophyletic clades")
public class CalibrationNode extends BEASTObject {
    public Node mrca;
    public CalibrationClade clade;
    public TaxonSet taxa;

    public CalibrationNode parent;
    public List<CalibrationNode> children = new ArrayList<>();

    public boolean isOverlapEdge; // true if overlaps parent interval
    public boolean isRoot;        // true if start of independent subtree

    // Fitted parameters
    public double mu;        // log-mean (lognormal)
    public double sigma2;    // log-variance (lognormal)

    double mEdge;
    double vEdge;                 // edge log-mean and log-variance increments

    public double alpha;     // Beta alpha
    public double beta;      // Beta beta

    // constructor used in CalibrationPrior
    public CalibrationNode(TreeInterface tree, CalibrationClade clade) {
        this.clade = clade;
        this.taxa = clade.getTaxa();
        this.mrca = getCommonAncestor(tree, clade.getTaxa());
    }

    // constructor used for CalibratedCPP
    public CalibrationNode(TreeInterface tree, TaxonSet taxa) {
        this.mrca = getCommonAncestor(tree, taxa);
        this.taxa = taxa;
    }

    public double getLower() {
        return clade != null ? clade.getLowerAge() : Double.NaN;
    }
    public double getUpper() {
        return clade != null ? clade.getUpperAge() : Double.NaN;
    }
    public double getCoverage() {
        return clade != null ? clade.getPCoverage() : Double.NaN;
    }

    public static CalibrationNode getByTaxa(List<CalibrationNode> nodes, TaxonSet taxa) {
        for (CalibrationNode node : nodes) {
            if (node.taxa != null && node.taxa.equals(taxa)) {
                return node;
            }
        }
        return null; // not found
    }

    private Node getCommonAncestor(TreeInterface tree, TaxonSet taxa) {
        List<Node> nodes = new ArrayList<>();
        for (Taxon t : taxa.getTaxonSet()) {
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
        String id = (taxa != null ? taxa.getID() : "noID");
        return String.format("CalibrationNode[%s, [%.3f, %.3f]]", id, getLower(), getUpper());
    }
    @Override
    public void initAndValidate() {

    }
}
