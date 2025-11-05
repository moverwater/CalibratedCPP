package calibratedcpp;

import calibratedcpp.model.BirthDeathModel;

import java.util.Arrays;
import java.util.List;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.alignment.Taxon;
import beast.base.inference.parameter.RealParameter;
import calibrationprior.CalibrationNode;
import org.junit.Test;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

public class CalibratedCoalescentPointProcessTest {

    private CalibratedCoalescentPointProcess cpp;
    private Tree tree;
    private TaxonSet cpABC;
    private TaxonSet cpDE;
    private TaxonSet cpABCDE;
    private TaxonSet cpHI;
    private TaxonSet cpHIJ;
    private TaxonSet cpFGHIJ;
    private BirthDeathModel birthDeath;
    private List<TaxonSet> calibrations;
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
        TaxonSet taxaABC = new TaxonSet(Arrays.asList(A, B, C));
        TaxonSet taxaDE = new TaxonSet(Arrays.asList(D, E));
        TaxonSet taxaABCDE = new TaxonSet(Arrays.asList(A, B, C, D, E));
        TaxonSet taxaHI = new TaxonSet(Arrays.asList(H, I));
        TaxonSet taxaHIJ = new TaxonSet(Arrays.asList(H, I, J));
        TaxonSet taxaFGHIJ = new TaxonSet(Arrays.asList(F, G, H, I, J));

        // Create CalibrationPoints
        cpABC = taxaABC;

        cpDE = taxaDE;

        cpABCDE = taxaABCDE;

        cpHI = taxaHI;

        cpHIJ = taxaHIJ;

        cpFGHIJ = taxaFGHIJ;

        calibrations = Arrays.asList(cpDE, cpABC, cpABCDE, cpHI, cpHIJ, cpFGHIJ);    }

    @Test
    public void calculateUnConditionedTreeLogLikelihood() {

        assertEquals(-25.05062, cpp.calculateUnConditionedTreeLogLikelihood(tree), 1e-4, "Unconditioned density of the tree is incorrect.");
        assertEquals(-25.05062 - Math.log(1 - Math.exp(birthDeath.calculateLogCDF(6.5))) + 2 * Math.log(1 - Math.exp(birthDeath.calculateLogCDF(6.0)))
                        - birthDeath.calculateLogDensity(6.0),
                rootConditionedCPP.calculateUnConditionedTreeLogLikelihood(tree), 1e-4, "Unconditioned density of the tree is incorrect."
                );
    }

    @Test
    public void calibrationForest(){
        List<CalibrationNode> calibrationGraph = cpp.buildCalibrationForest(tree, calibrations);

        for (CalibrationNode node : calibrationGraph) {
            if (!node.isRoot) {
                assertTrue( node.parent.taxa.getTaxonSet().size() > node.taxa.getTaxonSet().size());
                assertTrue(node.parent.mrca.getHeight() > node.mrca.getHeight());
            }
        }

        int n_roots = 0;
        for (CalibrationNode node : calibrationGraph) {
            n_roots += node.isRoot ? 1 : 0;
        }

        assertEquals(2, n_roots, "Calibration forest should have 2 roots but was " + n_roots);
    }

    @Test
    public void calculateLogMarginalDensityOfCalibrations() {
        List<CalibrationNode> calibrationGraph = cpp.buildCalibrationForest(tree, calibrations);

        assertEquals(birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                        Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
                        birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
                        Math.log(2.0) + // number of ways to arrange the clades
                        2 * Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.0))), // two terminating heights for the two sub-trees
                rootConditionedCPP.calculateLogMarginalDensityOfCalibrations(tree, calibrationGraph),1e-4, "Unconditioned density of the tree is incorrect.");

        assertEquals(birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                        Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
                        birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
                        Math.log(2.0) + // number of ways to arrange the clades
                        Math.log(Math.exp(birthDeath.calculateLogCDF(6.5)) - Math.exp(birthDeath.calculateLogCDF(5.0))) + Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5))), // distribution of the free node age and the terminating node age > origin
                cpp.calculateLogMarginalDensityOfCalibrations(tree, calibrationGraph), "MarginalDensity of the calibrations is incorrect.");

    }

//    @Test
//    public void calculateLogDensityOfSingleCalibration() {
//        assertEquals(birthDeath.calculateLogDensity(0.5), cpp.calculateLogDensityOfSingleCalibration(tree, cpHI, cpp.calibrationGraph), 1e-6, "Density for calibration HI is incorrect.");
//        assertEquals(birthDeath.calculateLogDensity(1.5), cpp.calculateLogDensityOfSingleCalibration(tree, cpDE, cpp.calibrationGraph), 1e-6, "Density for calibration DE is incorrect.");
//        assertEquals(birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0),
//                cpp.calculateLogDensityOfSingleCalibration(tree, cpABC, cpp.calibrationGraph), 1e-6, "Density for calibration ABC is incorrect.");
//        assertEquals(birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogCDF(2.5) + Math.log(2.0),
//                cpp.calculateLogDensityOfSingleCalibration(tree, cpHIJ, cpp.calibrationGraph), 1e-6, "Density for calibration HIJ is incorrect.");
//        assertEquals(3 * birthDeath.calculateLogCDF(4.0) + birthDeath.calculateLogDensity(4.0) + Math.log(4.0),
//                cpp.calculateLogDensityOfSingleCalibration(tree, cpABCDE, cpp.calibrationGraph), 1e-6, "Density for calibration ABCDE is incorrect.");
//
//        List<CalibrationNode> calibrationGraph = cpp.buildCalibrationForest(tree, calibrations);
//
//        List<CalibrationPoint> children = calibrationGraph.getOrDefault(cpABCDE, new ArrayList<>());
//        assertEquals(cpABC, children.get(0), "Index of calibration ABC is incorrect.");
//
//        assertEquals(birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0),
//                cpp.calculateLogDensityOfSingleCalibration(tree, cpHIJ, calibrationGraph), 1e-6, "Density for calibration HIJ is incorrect.");
//        assertEquals(birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0),
//                cpp.calculateLogDensityOfSingleCalibration(tree, cpABCDE, calibrationGraph), 1e-6, "Density for calibration ABCDE is incorrect.");
//        assertEquals(birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
//                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))),
//                cpp.calculateLogDensityOfSingleCalibration(tree, cpFGHIJ, calibrationGraph), 1e-6, "Density for calibration FGHIJ is incorrect.");
//    }
//
//    @Test
//    public void calculateTreeLogLikelihood() {
//        cpp = new CalibratedCoalescentPointProcess();
//        cpp.initByName("tree", tree,
//                "origin", new RealParameter("6.5"),
//                "calibrations", calibrations,
//                "treeModel", birthDeath);
//
//        rootConditionedCPP = new CalibratedCoalescentPointProcess();
//        rootConditionedCPP.initByName("tree", tree,
//                "origin", new RealParameter("6.5"),
//                "calibrations", calibrations,
//                "conditionOnRoot", true,
//                "conditionOnCalibrations", true,
//                "treeModel", birthDeath);
//
//        assertEquals(-25.05062 - // unconditioned tree density
//                        (birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
//                                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
//                                birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
//                                Math.log(2.0) + // number of ways to arrange the clades
//                                Math.log(Math.exp(birthDeath.calculateLogCDF(cpp.origin)) - Math.exp(birthDeath.calculateLogCDF(5.0))) + Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5)))), // distribution of the free node age and the terminating node age > origin,
//                cpp.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");
//
//        assertEquals(-25.05062 - birthDeath.calculateLogDensity(6.0) - Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5))) - // unconditioned tree density
//                        (birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
//                                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
//                                birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
//                                Math.log(2.0)), // number of ways to arrange the clades
//                rootConditionedCPP.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");
//    }
//
}