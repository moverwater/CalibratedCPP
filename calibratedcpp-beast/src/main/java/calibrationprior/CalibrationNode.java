package calibrationprior;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

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
    public CalibrationNode(Node mrca, CalibrationClade clade) {
        this.mrca = mrca;
        this.clade = clade;
        this.taxa = clade.getTaxa();
    }

    // constructor used for CalibratedCPP
    public CalibrationNode(Node mrca, TaxonSet taxa) {
        this.mrca = mrca;
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

    @Override
    public String toString() {
        String id = (taxa != null ? taxa.getID() : "noID");
        return String.format("CalibrationNode[%s, [%.3f, %.3f]]", id, getLower(), getUpper());
    }
    @Override
    public void initAndValidate() {

    }
}
