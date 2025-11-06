package calibrationprior;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;

import java.util.*;

/**
 * @author Marcus Overwater
 */

@Description("A forest of calibration nodes.")
public class CalibrationForest extends BEASTObject {

    private final List<CalibrationNode> allNodes = new ArrayList<>();
    private final List<CalibrationNode> roots = new ArrayList<>();

    // --- Private constructor (used only internally) ---
    private CalibrationForest() {}

    // --- Static factory for clade-based forest ---
    public static CalibrationForest buildFromClades(TreeInterface tree, List<CalibrationClade> clades) {
        CalibrationForest forest = new CalibrationForest();
        for (CalibrationClade clade : clades) {
            forest.allNodes.add(new CalibrationNode(tree, clade));
        }
        forest.buildInclusionForest();
        return forest;
    }

    // --- Static factory for taxa-based forest (no bounds) ---
    public static CalibrationForest buildFromTaxonSets(TreeInterface tree, List<TaxonSet> taxaSets) {
        CalibrationForest forest = new CalibrationForest();
        for (TaxonSet tset : taxaSets) {
            forest.allNodes.add(new CalibrationNode(tree, tset));
        }
        forest.buildInclusionForest();
        return forest;
    }

    // --- Core forest structure builder ---
    private void buildInclusionForest() {
        roots.clear(); // clear old roots

        // Precompute taxa sets for each node
        Map<CalibrationNode, Set<String>> nodeTaxaMap = new HashMap<>();
        for (CalibrationNode node : allNodes) {
            Set<String> taxaIds = new HashSet<>();
            for (var t : node.taxa.getTaxonSet()) {
                taxaIds.add(t.getID());
            }
            nodeTaxaMap.put(node, taxaIds);
        }

        // Validate nested or disjoint
        for (int i = 0; i < allNodes.size(); i++) {
            CalibrationNode a = allNodes.get(i);
            Set<String> aTaxa = nodeTaxaMap.get(a);
            for (int j = i + 1; j < allNodes.size(); j++) {
                CalibrationNode b = allNodes.get(j);
                Set<String> bTaxa = nodeTaxaMap.get(b);

                boolean aContainsB = aTaxa.containsAll(bTaxa);
                boolean bContainsA = bTaxa.containsAll(aTaxa);

                if (!aContainsB && !bContainsA && !Collections.disjoint(aTaxa, bTaxa)) {
                    throw new IllegalArgumentException(
                            "Calibration clades must be nested or disjoint: " +
                                    a.taxa.getID() + " and " + b.taxa.getID()
                    );
                }
            }
        }

        // Build parent-child relations
        for (CalibrationNode child : allNodes) {
            CalibrationNode bestParent = null;
            Set<String> childTaxa = nodeTaxaMap.get(child);

            for (CalibrationNode candidate : allNodes) {
                if (candidate == child) continue;
                Set<String> candidateTaxa = nodeTaxaMap.get(candidate);

                if (candidateTaxa.containsAll(childTaxa)) {
                    // pick the most specific parent (smallest superset)
                    if (bestParent == null || nodeTaxaMap.get(bestParent).containsAll(candidateTaxa)) {
                        bestParent = candidate;
                    }
                }
            }

            child.parent = bestParent;
            child.isRoot = (bestParent == null);

            if (bestParent != null) {
                bestParent.children.add(child);
            } else {
                roots.add(child);
            }
        }
    }

    // --- Lookup methods ---
    public CalibrationNode getCalibrationNodeFromTaxonSet(TaxonSet taxonSet) {
        for (CalibrationNode node : allNodes) {
            if (node.taxa != null && node.taxa.getTaxonSet().equals(taxonSet.getTaxonSet())) {
                return node;
            }
        }
        return null;
    }

    // --- Accessors ---
    public List<CalibrationNode> getAllNodes() { return allNodes; }
    public List<CalibrationNode> getRoots() { return roots; }

    @Override
    public void initAndValidate() {
        buildInclusionForest();
    }
}
