package calibration;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import calibrationprior.CalibrationCladePrior;

import java.util.*;

/**
 * A node in a calibration forest.
 *
 * @author Marcus Overwater
 */

@Description("A calibration node is a node in a calibration forest.")
public class CalibrationNode extends BEASTObject {
    public TaxonSet taxa;

    /** Non-null when this node was created from a {@link CalibrationCladePrior}. */
    public CalibrationCladePrior prior;

    public CalibrationNode parent;
    public List<CalibrationNode> children = new ArrayList<>();

    public boolean isRoot;        // true if start of independent subtree

    // Cached leaf map: taxon name -> node number (stable since nodes are stored in parallel arrays)
    private static final WeakHashMap<TreeInterface, Map<String, Integer>> leafMapCache = new WeakHashMap<>();

    public CalibrationNode(TaxonSet taxa) {
        this.taxa = taxa;
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

    public CalibrationCladePrior getCalibrationCladePrior() {
        if (prior == null) {
            throw new IllegalStateException("CalibrationNode does not hold a CalibrationCladePrior");
        }
        return prior;
    }

    // --- Utility ---
    public Node getCommonAncestor(TreeInterface tree) {
        // Cache node numbers per tree - O(1) lookup vs BEAST's O(n) getTaxonIndex
        Map<String, Integer> leafMap = leafMapCache.get(tree);
        if (leafMap == null) {
            leafMap = new HashMap<>();
            for (Node n : tree.getExternalNodes()) {
                leafMap.put(n.getID(), n.getNr());
            }
            leafMapCache.put(tree, leafMap);
        }

        List<Node> nodes = new ArrayList<>();
        for (Taxon t : this.taxa.getTaxonSet()) {
            Integer nr = leafMap.get(t.getID());
            if (nr != null) {
                nodes.add(tree.getNode(nr));
            }
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
        while (a != null) {
            ancA.add(a);
            a = a.getParent();
        }
        while (b != null) {
            if (ancA.contains(b)) return b;
            b = b.getParent();
        }
        return null;
    }

    @Override
    public String toString() {
        return String.format("CalibrationNode[%s]", taxa != null ? taxa.getID() : "noID");
    }

    @Override
    public void initAndValidate() {}
}
