package calibration;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A forest of CalibrationNodes
 *
 * @author Marcus Overwater
 */

@Description("A forest of calibration nodes.")
public class CalibrationForest extends BEASTObject {

    private final List<CalibrationNode> allNodes = new ArrayList<>();
    private final List<CalibrationNode> roots = new ArrayList<>();

    // --- Constructor ---
    public CalibrationForest(List<? extends CalibrationClade> clades) {
        for (CalibrationClade clade : clades) {
            allNodes.add(new CalibrationNode(clade));
        }
        buildInclusionForest();
    }

    // --- Core forest structure builder ---
    private void buildInclusionForest() {
        roots.clear();

        // Precompute taxa sets
        Map<CalibrationNode, Set<String>> nodeTaxaMap = new HashMap<>();
        for (CalibrationNode node : allNodes) {
            Set<String> taxaIds = node.taxa.getTaxonSet().stream()
                    .map(Taxon::getID)
                    .collect(Collectors.toSet());
            nodeTaxaMap.put(node, taxaIds);
        }

        // --- Validate nesting/disjointness ---
        for (int i = 0; i < allNodes.size(); i++) {
            CalibrationNode a = allNodes.get(i);
            Set<String> aTaxa = nodeTaxaMap.get(a);

            for (int j = i + 1; j < allNodes.size(); j++) {
                CalibrationNode b = allNodes.get(j);
                Set<String> bTaxa = nodeTaxaMap.get(b);

                boolean aContainsB = aTaxa.containsAll(bTaxa);
                boolean bContainsA = bTaxa.containsAll(aTaxa);
                boolean disjoint = Collections.disjoint(aTaxa, bTaxa);

                // identical taxon sets not allowed
                if (aTaxa.equals(bTaxa)) {
                    throw new IllegalArgumentException("Duplicate calibration clade: " + a.taxa.getID());
                }

                if (!aContainsB && !bContainsA && !disjoint) {
                    throw new IllegalArgumentException(
                            "Calibration clades must be nested or disjoint: " +
                                    a.taxa.getID() + " and " + b.taxa.getID());
                }
            }
        }

        // --- Build parent-child relations ---
        for (CalibrationNode child : allNodes) {
            CalibrationNode bestParent = null;
            Set<String> childTaxa = nodeTaxaMap.get(child);

            for (CalibrationNode candidate : allNodes) {
                if (candidate == child) continue;

                Set<String> candidateTaxa = nodeTaxaMap.get(candidate);
                if (candidateTaxa.containsAll(childTaxa)) {
                    // most specific superset (smallest one)
                    if (bestParent == null ||
                            (nodeTaxaMap.get(bestParent).containsAll(candidateTaxa)
                                    && !nodeTaxaMap.get(candidate).equals(childTaxa))) {
                        bestParent = candidate;
                    }
                }
            }

            // Prevent self or circular assignment
            if (bestParent == child) {
                throw new IllegalStateException("Self-parent detected for " + child);
            }

            child.parent = bestParent;
            child.isRoot = (bestParent == null);

            if (bestParent != null) {
                bestParent.children.add(child);
            } else {
                roots.add(child);
            }
        }

        // Final sanity check: detect cycles
        for (CalibrationNode node : allNodes) {
            checkForCycle(node, new HashSet<>());
        }
    }

    // Helper to detect cycles
    private void checkForCycle(CalibrationNode node, Set<CalibrationNode> visited) {
        if (node == null) return;
        if (!visited.add(node)) {
            throw new IllegalStateException("Cycle detected involving: " + node);
        }
        if (node.parent != null) {
            checkForCycle(node.parent, visited);
        }
    }

    // --- Lookup methods ---
    public CalibrationNode getCalibrationNodeFromTaxonSet(TaxonSet taxonSet) {
        for (CalibrationNode node : allNodes) {
            Set<String> nodeIds = node.taxa.getTaxonSet().stream().map(Taxon::getID).collect(Collectors.toSet());
            Set<String> targetIds = taxonSet.getTaxonSet().stream().map(Taxon::getID).collect(Collectors.toSet());
            if (nodeIds.equals(targetIds)) return node;
        }
        return null;
    }

    public CalibrationNode getCalibrationNodeFromCalibrationClade(CalibrationClade calibrationClade) {
        return getCalibrationNodeFromTaxonSet(calibrationClade.getTaxa());
    }

    // --- Accessors ---
    public List<CalibrationNode> getAllNodes() { return allNodes; }
    public List<CalibrationNode> getRoots() { return roots; }

    @Override
    public void initAndValidate() {
        buildInclusionForest();
    }
}
