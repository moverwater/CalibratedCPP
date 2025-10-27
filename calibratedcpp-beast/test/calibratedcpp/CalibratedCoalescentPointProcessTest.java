package calibratedcpp;

import calibratedcpp.model.BirthDeathModel;

import java.util.Arrays;
import java.util.List;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.alignment.Taxon;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;

import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

public class CalibratedCoalescentPointProcessTest {

    private CalibratedCoalescentPointProcess cpp;
    private Tree tree;
    private CalibrationPoint cpABC;
    private CalibrationPoint cpDE;
    private CalibrationPoint cpABCDE;
    private CalibrationPoint cpHI;
    private CalibrationPoint cpHIJ;
    private CalibrationPoint cpFGHIJ;
    private BirthDeathModel birthDeath;
    private List<CalibrationPoint> calibrations;
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

        // Create calibration distributions
        ParametricDistribution distABC = new Uniform();
        distABC.initByName("lower", 2.9, "upper", 3.1);

        ParametricDistribution distDE = new Uniform();
        distDE.initByName("lower", 1.4, "upper", 1.6);

        ParametricDistribution distABCDE = new Uniform();
        distABCDE.initByName("lower", 3.9, "upper", 4.1);

        ParametricDistribution distHI = new Uniform();
        distHI.initByName("lower", 0.4, "upper", 0.6);

        ParametricDistribution distHIJ = new Uniform();
        distHIJ.initByName("lower", 2.4, "upper", 2.6);

        ParametricDistribution distFGHIJ = new Uniform();
        distFGHIJ.initByName("lower", 4.9, "upper", 5.1);

        // Create CalibrationPoints
        cpABC = new CalibrationPoint();
        cpABC.initByName("taxonset", taxaABC, "distr", distABC);

        cpDE = new CalibrationPoint();
        cpDE.initByName("taxonset", taxaDE, "distr", distDE);

        cpABCDE = new CalibrationPoint();
        cpABCDE.initByName("taxonset", taxaABCDE, "distr", distABCDE);

        cpHI = new CalibrationPoint();
        cpHI.initByName("taxonset", taxaHI, "distr", distHI);

        cpHIJ = new CalibrationPoint();
        cpHIJ.initByName("taxonset", taxaHIJ, "distr", distHI);

        cpFGHIJ = new CalibrationPoint();
        cpFGHIJ.initByName("taxonset", taxaFGHIJ, "distr", distHI);

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
    public void calculateLogMarginalDensityOfCalibrations() {
        Map<CalibrationPoint, List<CalibrationPoint>> calibrationGraph = cpp.buildNestingDAG(calibrations, tree);

        assertEquals(birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                        Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
                        birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
                        Math.log(2.0) + // number of ways to arrange the clades
                        Math.log(Math.exp(birthDeath.calculateLogCDF(cpp.origin)) - Math.exp(birthDeath.calculateLogCDF(5.0))) + Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5))), // distribution of the free node age and the terminating node age > origin
                cpp.calculateLogMarginalDensityOfCalibrations(tree, calibrations, calibrationGraph));

        assertEquals(birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
                birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
                Math.log(2.0) + // number of ways to arrange the clades
                2 * Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.0))), // two terminating heights for the two sub-trees
                rootConditionedCPP.calculateLogMarginalDensityOfCalibrations(tree, calibrations, calibrationGraph),1e-4, "Unconditioned density of the tree is incorrect.");
    }

    @Test
    public void calculateLogDensityOfSingleCalibration() {
        assertEquals(birthDeath.calculateLogDensity(0.5), cpp.calculateLogDensityOfSingleCalibration(tree, cpHI, cpp.calibrationGraph), 1e-6, "Density for calibration HI is incorrect.");
        assertEquals(birthDeath.calculateLogDensity(1.5), cpp.calculateLogDensityOfSingleCalibration(tree, cpDE, cpp.calibrationGraph), 1e-6, "Density for calibration DE is incorrect.");
        assertEquals(birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0),
                cpp.calculateLogDensityOfSingleCalibration(tree, cpABC, cpp.calibrationGraph), 1e-6, "Density for calibration ABC is incorrect.");
        assertEquals(birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogCDF(2.5) + Math.log(2.0),
                cpp.calculateLogDensityOfSingleCalibration(tree, cpHIJ, cpp.calibrationGraph), 1e-6, "Density for calibration HIJ is incorrect.");
        assertEquals(3 * birthDeath.calculateLogCDF(4.0) + birthDeath.calculateLogDensity(4.0) + Math.log(4.0),
                cpp.calculateLogDensityOfSingleCalibration(tree, cpABCDE, cpp.calibrationGraph), 1e-6, "Density for calibration ABCDE is incorrect.");

        Map<CalibrationPoint, List<CalibrationPoint>> calibrationGraph = cpp.buildNestingDAG(calibrations, tree);

        List<CalibrationPoint> children = calibrationGraph.getOrDefault(cpABCDE, new ArrayList<>());
        assertEquals(cpABC, children.get(0), "Index of calibration ABC is incorrect.");

        assertEquals(birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0),
                cpp.calculateLogDensityOfSingleCalibration(tree, cpHIJ, calibrationGraph), 1e-6, "Density for calibration HIJ is incorrect.");
        assertEquals(birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0),
                cpp.calculateLogDensityOfSingleCalibration(tree, cpABCDE, calibrationGraph), 1e-6, "Density for calibration ABCDE is incorrect.");
        assertEquals(birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))),
                cpp.calculateLogDensityOfSingleCalibration(tree, cpFGHIJ, calibrationGraph), 1e-6, "Density for calibration FGHIJ is incorrect.");
    }

    @Test
    public void calculateTreeLogLikelihood() {
        cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "origin", new RealParameter("6.5"),
                "calibrations", calibrations,
                "treeModel", birthDeath);

        rootConditionedCPP = new CalibratedCoalescentPointProcess();
        rootConditionedCPP.initByName("tree", tree,
                "origin", new RealParameter("6.5"),
                "calibrations", calibrations,
                "conditionOnRoot", true,
                "conditionOnCalibrations", true,
                "treeModel", birthDeath);

        assertEquals(-25.05062 - // unconditioned tree density
                        (birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
                                birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
                                Math.log(2.0) + // number of ways to arrange the clades
                                Math.log(Math.exp(birthDeath.calculateLogCDF(cpp.origin)) - Math.exp(birthDeath.calculateLogCDF(5.0))) + Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5)))), // distribution of the free node age and the terminating node age > origin,
                cpp.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");

        assertEquals(-25.05062 - birthDeath.calculateLogDensity(6.0) - Math.log1p(-Math.exp(birthDeath.calculateLogCDF(6.5))) - // unconditioned tree density
                        (birthDeath.calculateLogDensity(5.0) + birthDeath.calculateLogDensity(2.5) + birthDeath.calculateLogDensity(0.5) + Math.log(2.0) +
                                Math.log(4 * (Math.exp(birthDeath.calculateLogCDF(5.0)) - Math.exp(birthDeath.calculateLogCDF(2.5))) + 2 * Math.exp(birthDeath.calculateLogCDF(5.0))) + // density of FGHIJ
                                birthDeath.calculateLogDensity(4.0) + birthDeath.calculateLogDensity(3.0) + birthDeath.calculateLogDensity(1.5) + birthDeath.calculateLogCDF(3.0) + Math.log(2.0) + Math.log(2.0) + // density of ABCDE
                                Math.log(2.0)), // number of ways to arrange the clades
                rootConditionedCPP.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");
    }

    @Test
    public void isNested() {
        assertTrue(cpp.isNested(cpABC, cpABCDE, tree), "A,B should be nested inside A,B,C");
        assertTrue(cpp.isNested(cpHI, cpFGHIJ, tree), "HI should be nested inside F,G,H,I,J");

        assertFalse(cpp.isNested(cpABC, cpHI, tree), "A,B,C should not be nested inside H,I");
        assertFalse(cpp.isNested(cpDE, cpABC, tree), "D,E should not be nested inside A,B,C");
        assertFalse(cpp.isNested(cpHIJ, cpHI, tree), "H,I,J should be nested inside H,I");
    }


    @Test
    public void buildNestingDAG() {
        List<CalibrationPoint> calibrations = Arrays.asList(cpABC, cpABCDE, cpDE, cpFGHIJ, cpHI, cpHIJ);
        Map<CalibrationPoint, List<CalibrationPoint>> graph = cpp.buildNestingDAG(calibrations, tree);

        assertTrue(graph.get(cpABCDE).contains(cpABC), "cpABC should have cpAB as child");
        assertTrue(graph.get(cpABCDE).contains(cpDE), "cpACDE should have cpDE as child");
        assertTrue(graph.get(cpFGHIJ).contains(cpHIJ), "cpFGHIJ should have cpHIJ as child");
        assertTrue(graph.get(cpHIJ).contains(cpHI), "cpHIJ should have cpHI as child");

        assertFalse(graph.get(cpFGHIJ).contains(cpHI), "cpFGHIJ should not have cpHI as child");
        assertFalse(graph.get(cpFGHIJ).contains(cpABC), "cpFGHIJ should not have cpABC as child");
        assertFalse(graph.get(cpHI).contains(cpHIJ), "cpHI should not have cpHIJ as child");

        assertTrue(graph.get(cpABC).isEmpty(), "cpABC should have no children");
        assertTrue(graph.get(cpHI).isEmpty(), "cpHI should have no children");
    }

    @Test
    public void postOrderTopologicalSort() {
        List<CalibrationPoint> calibrations = Arrays.asList(cpABC, cpABCDE, cpDE, cpFGHIJ, cpHI, cpHIJ);
        List<CalibrationPoint> sorted = cpp.postOrderTopologicalSort(tree, calibrations);

        int indexABC = sorted.indexOf(cpABC);
        int indexABCDE = sorted.indexOf(cpABCDE);
        int indexHIJ = sorted.indexOf(cpHIJ);
        int indexHI = sorted.indexOf(cpHI);
        int indexFGHIJ = sorted.indexOf(cpFGHIJ);
        int indexDE = sorted.indexOf(cpDE);

        // Children must come before parent
        assertTrue(indexABC < indexABCDE, "cpAB should come before cpABC");
        assertTrue(indexDE < indexABCDE, "cDE should come before cpABCDE");
        assertTrue(indexDE > indexABC, "cDE should come before cpABC");
        assertTrue(indexHI < indexHIJ, "cpAC should come before cpABC");
        assertTrue(indexFGHIJ > indexABCDE, "cpFGHIJ should come before cpABCDE");
    }
}