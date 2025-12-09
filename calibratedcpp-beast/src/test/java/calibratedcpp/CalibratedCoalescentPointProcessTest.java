package calibratedcpp;

import calibratedcpp.model.BirthDeathModel;

import java.util.Arrays;
import java.util.List;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.alignment.Taxon;
import beast.base.inference.parameter.RealParameter;
import calibration.CalibrationClade;
import calibration.CalibrationForest;
import calibration.CalibrationNode;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class CalibratedCoalescentPointProcessTest {

    private CalibratedCoalescentPointProcess cpp;
    private final Tree tree;
    TaxonSet taxaABC;
    TaxonSet taxaDE;
    TaxonSet taxaABCDE;
    TaxonSet taxaHI;
    TaxonSet taxaHIJ;
    TaxonSet taxaFGHIJ;
    private final BirthDeathModel birthDeath;
    private final List<CalibrationClade> calibrations;
    private CalibratedCoalescentPointProcess rootConditionedCPP;

    public CalibratedCoalescentPointProcessTest() {
        tree = new TreeParser();
        tree.initByName("newick", "((((A:2,B:2):1,C:3):1,(D:1.5,E:1.5):2.5):2,((F:1,G:1):4,((H:0.5,I:0.5):2,J:2.5):2.5):1):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);

        birthDeath = new BirthDeathModel();
        birthDeath.initByName("birthRate", new RealParameter("3.0"),
                "deathRate", new RealParameter("2.0"),
                "rho", new RealParameter("0.1")
        );

        cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "treeModel", birthDeath,
                "origin", new RealParameter("6.5")
        );

        rootConditionedCPP = new CalibratedCoalescentPointProcess();
        rootConditionedCPP.initByName("tree", tree,
                "treeModel", birthDeath,
                "conditionOnRoot", true
        );

        // Create taxa objects
        Taxon A = new Taxon();
        A.setID("A");
        Taxon B = new Taxon();
        B.setID("B");
        Taxon C = new Taxon();
        C.setID("C");
        Taxon D = new Taxon();
        D.setID("D");
        Taxon E = new Taxon();
        E.setID("E");
        Taxon F = new Taxon();
        F.setID("F");
        Taxon G = new Taxon();
        G.setID("G");
        Taxon H = new Taxon();
        H.setID("H");
        Taxon I = new Taxon();
        I.setID("I");
        Taxon J = new Taxon();
        J.setID("J");

        // Create TaxonSets
        taxaABC = new TaxonSet(Arrays.asList(A, B, C));
        taxaDE = new TaxonSet(Arrays.asList(D, E));
        taxaABCDE = new TaxonSet(Arrays.asList(A, B, C, D, E));
        taxaHI = new TaxonSet(Arrays.asList(H, I));
        taxaHIJ = new TaxonSet(Arrays.asList(H, I, J));
        taxaFGHIJ = new TaxonSet(Arrays.asList(F, G, H, I, J));

        CalibrationClade calibrationABC = new CalibrationClade();
        calibrationABC.initByName("taxa",taxaABC);
        CalibrationClade calibrationDE = new CalibrationClade();
        calibrationDE.initByName("taxa",taxaDE);
        CalibrationClade calibrationABCDE = new CalibrationClade();
        calibrationABCDE.initByName("taxa",taxaABCDE);
        CalibrationClade calibrationHI = new CalibrationClade();
        calibrationHI.initByName("taxa",taxaHI);
        CalibrationClade calibrationHIJ = new CalibrationClade();
        calibrationHIJ.initByName("taxa",taxaHIJ);
        CalibrationClade calibrationFGHIJ = new CalibrationClade();
        calibrationFGHIJ.initByName("taxa",taxaFGHIJ);

        calibrations = Arrays.asList(calibrationDE, calibrationABC, calibrationABCDE, calibrationHI, calibrationHIJ, calibrationFGHIJ);    }

    @Test
    public void calculateUnConditionedTreeLogLikelihood() {

        int n = tree.getLeafNodeCount();

        double log_n_factorial = 0;
        for (int i = 1 ; i <= n; i++){
            log_n_factorial += Math.log(i);
        }

        assertEquals(-25.05062 + (n - 1) * Math.log(2) - log_n_factorial,
                cpp.calculateUnConditionedTreeLogLikelihood(tree), 1e-4, "Unconditioned density of the tree is incorrect.");
    }

    @Test
    public void calibrationForest(){
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        List<CalibrationNode> calibrationNodes = calibrationForest.getAllNodes();

        for (CalibrationNode node : calibrationNodes) {
            if (!node.isRoot) {
                assertTrue( node.parent.getCalibrationClade().getTaxa().getTaxonCount() > node.getCalibrationClade().getTaxa().getTaxonCount());
                assertTrue(node.parent.getCommonAncestor(tree).getHeight() > node.getCommonAncestor(tree).getHeight());
            }
        }

        int n_roots = 0;
        for (CalibrationNode node : calibrationNodes) {
            n_roots += node.isRoot ? 1 : 0;
        }

        assertEquals(2, n_roots, "Calibration forest should have 2 roots but was " + n_roots);
    }

    @Test
    public void calculateLogMarginalDensityOfCalibrations() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        List<CalibrationNode> calibrationNodes = calibrationForest.getAllNodes();

        CalibrationNode cpFGHIJ = CalibrationNode.getByTaxa(calibrationNodes, taxaFGHIJ);
        CalibrationNode cpABCDE = CalibrationNode.getByTaxa(calibrationNodes, taxaABCDE);

        double labellings = 0;
        for (int i = 1; i <= 5 ; i++ ){
            labellings += Math.log(i);
        }
        labellings = 2 * labellings;
        for (int i = 1; i <= 10 ; i++ ){
            labellings -= Math.log(i);
        }

        assert cpABCDE != null;
        assert cpFGHIJ != null;
        assertEquals(rootConditionedCPP.computeCalibrationDensity(tree, cpABCDE) +
                        rootConditionedCPP.computeCalibrationDensity(tree, cpFGHIJ) +
                        Math.log(2.0) + birthDeath.calculateLogDensity(6.0) + labellings,
                rootConditionedCPP.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest),1e-4,
                "Marginal density of the calibrations and root is incorrect.");

        assertEquals(cpp.computeCalibrationDensity(tree, cpABCDE) + cpp.computeCalibrationDensity(tree, cpFGHIJ) +
                Math.log(2.0) + birthDeath.calculateLogCDF(6.5) + Math.log1p(-Math.exp(birthDeath.calculateLogCDF(5.0) - birthDeath.calculateLogCDF(6.5))) +
                labellings,
                cpp.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest), 1e-4,
                "MarginalDensity of the calibrations is incorrect.");
    }

    @Test
    public void computeCalibrationDensity() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        List<CalibrationNode> calibrationNodes = calibrationForest.getAllNodes();

        CalibrationNode cpHI = CalibrationNode.getByTaxa(calibrationNodes, taxaHI);
        CalibrationNode cpDE = CalibrationNode.getByTaxa(calibrationNodes, taxaDE);
        CalibrationNode cpABC = CalibrationNode.getByTaxa(calibrationNodes, taxaABC);
        CalibrationNode cpHIJ = CalibrationNode.getByTaxa(calibrationNodes, taxaHIJ);
        CalibrationNode cpFGHIJ = CalibrationNode.getByTaxa(calibrationNodes, taxaFGHIJ);
        CalibrationNode cpABCDE = CalibrationNode.getByTaxa(calibrationNodes, taxaABCDE);

        assert cpHI != null;
        assertEquals(birthDeath.calculateLogDensity(0.5),
                cpp.computeCalibrationDensity(tree, cpHI), 1e-6, "Density for calibration HI is incorrect.");
        assert cpDE != null;
        assertEquals(birthDeath.calculateLogDensity(1.5),
                cpp.computeCalibrationDensity(tree, cpDE), 1e-6, "Density for calibration DE is incorrect.");
        assert cpABC != null;
        assertEquals(birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0),
                cpp.computeCalibrationDensity(tree, cpABC), 1e-6, "Density for calibration ABC is incorrect.");
        assert cpHIJ != null;
        assertEquals(birthDeath.calculateLogDensity(0.5) + birthDeath.calculateLogDensity(2.5) + Math.log(2.0)
                         - (Math.log(3.0)),
                cpp.computeCalibrationDensity(tree, cpHIJ), 1e-6, "Density for calibration HIJ is incorrect.");
        assert cpABCDE != null;
        assertEquals(birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) +
                        birthDeath.calculateLogDensity(1.5) +
                        birthDeath.calculateLogDensity(4.0) +
                        2 * Math.log(2.0) - (Math.log(5.0) + Math.log(4.0)),
                cpp.computeCalibrationDensity(tree, cpABCDE), 1e-6, "Density for calibration ABCDE is incorrect.");
        assert cpFGHIJ != null;
        assertEquals(birthDeath.calculateLogDensity(0.5) + birthDeath.calculateLogDensity(2.5) + Math.log(2.0)
                        - (Math.log(3.0)) + // density of HIJ
                birthDeath.calculateLogDensity(5.0) +
                Math.log(4.0 * (Math.exp(birthDeath.calculateLogCDF(5.0))-Math.exp(birthDeath.calculateLogCDF(2.5)))
                + 2.0 * Math.exp(birthDeath.calculateLogCDF(5.0))) + Math.log(2.0) - (Math.log(5.0) + Math.log(4.0)),
                cpp.computeCalibrationDensity(tree, cpFGHIJ), 1e-6, "Density for calibration FGHIJ is incorrect.");
    }

    @Test
    public void calculateTreeLogLikelihood() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "origin", new RealParameter("6.5"),
                "calibrations", calibrations,
                "treeModel", birthDeath);

        assertEquals(cpp.calculateUnConditionedTreeLogLikelihood(tree) -
                        cpp.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest) -
                Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5))),
                cpp.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");

        rootConditionedCPP = new CalibratedCoalescentPointProcess();
        rootConditionedCPP.initByName("tree", tree,
                "origin", new RealParameter("6.5"),
                "calibrations", calibrations,
                "conditionOnRoot", true,
                "conditionOnCalibrations", true,
                "treeModel", birthDeath);

        assertEquals(rootConditionedCPP.calculateUnConditionedTreeLogLikelihood(tree) -
                rootConditionedCPP.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest) -
                Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.0))),
                rootConditionedCPP.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");

        CalibratedCoalescentPointProcess nonmonophyleticCPP = new CalibratedCoalescentPointProcess();
        CalibrationClade badClade = new CalibrationClade();
        badClade.initByName("taxa", new TaxonSet(Arrays.asList(
                new Taxon("A"), new Taxon("B"), new Taxon("J")
        )));
        nonmonophyleticCPP.initByName("tree", tree,
                "origin", new RealParameter("6.5"),
                "calibrations", badClade,
                "conditionOnRoot", true,
                "conditionOnCalibrations", true,
                "treeModel", birthDeath);

        assertEquals(Double.NEGATIVE_INFINITY, nonmonophyleticCPP.calculateTreeLogLikelihood(tree), 1e-4);
    }

}