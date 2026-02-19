package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class CalibratedBirthDeathModelTest {

    private final CalibratedBirthDeathModel birthDeathModelSuperCrit;
    private final CalibratedBirthDeathModel birthDeathModelSubCrit;
    private final CalibratedBirthDeathModel birthDeathModelCrit;

    public CalibratedBirthDeathModelTest() {
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:1,B:1):1,C:2);",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        birthDeathModelSuperCrit = new CalibratedBirthDeathModel();
        birthDeathModelSuperCrit.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealParameter("2.0"),
                "deathRate", new RealParameter("1.0"),
                "rho", new RealParameter("0.1"));

        birthDeathModelSubCrit = new CalibratedBirthDeathModel();
        birthDeathModelSubCrit.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealParameter("1.0"),
                "deathRate", new RealParameter("2.0"),
                "rho", new RealParameter("0.1"));

        birthDeathModelCrit = new CalibratedBirthDeathModel();
        birthDeathModelCrit.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealParameter("1.0"),
                "diversificationRate", new RealParameter("0.0"),
                "rho", new RealParameter("1.0"));
    }

    @Test
    public void testBirthDeathParameters() {
        assertNotNull(birthDeathModelSuperCrit.birthRateInput.get(), "Birth rate should be initialized");
        assertNotNull(birthDeathModelSuperCrit.deathRateInput.get(), "Death rate should be initialized");
        assertEquals(2.0, birthDeathModelSuperCrit.birthRate, 1e-6, "Birth rate value incorrect");
        assertEquals(1.0, birthDeathModelSuperCrit.deathRate, 1e-6, "Death rate value incorrect");
        assertTrue(birthDeathModelSuperCrit.diversificationRate > 0);

        assertNotNull(birthDeathModelSubCrit.birthRateInput.get(), "Birth rate should be initialized");
        assertNotNull(birthDeathModelSubCrit.deathRateInput.get(), "Death rate should be initialized");
        assertEquals(1.0, birthDeathModelSubCrit.birthRate, 1e-6, "Birth rate value incorrect");
        assertEquals(2.0, birthDeathModelSubCrit.deathRate, 1e-6, "Death rate value incorrect");
        assertTrue(birthDeathModelSubCrit.diversificationRate < 0);

        assertNotNull(birthDeathModelCrit.birthRateInput.get(), "Birth rate should be initialized");
        assertNotNull(birthDeathModelCrit.diversificationRateInput.get(), "Diversification rate should be initialized");
        assertEquals(1.0, birthDeathModelCrit.birthRate, 1e-6, "Birth rate value incorrect");
        assertEquals(1.0, birthDeathModelCrit.deathRate, 1e-6, "Death rate value incorrect");
        assertEquals(0.0, birthDeathModelCrit.diversificationRate);
    }
    
    @Test
    public void calculateLogeNodeAgeDensityTest() {
        assertEquals(-8.390925, birthDeathModelSuperCrit.calculateLogNodeAgeDensity(10), 1e-4,
                "Log density incorrect.");
        assertEquals(-3.443752, birthDeathModelSuperCrit.calculateLogNodeAgeDensity(5), 1e-4,
                "Log density incorrect.");
        assertEquals(-1.200227, birthDeathModelSuperCrit.calculateLogNodeAgeDensity(1), 1e-4,
                "Log density incorrect.");

        assertEquals(-12.4932, birthDeathModelSubCrit.calculateLogNodeAgeDensity(10), 1e-4,
                "Log density incorrect.");
        assertEquals(-7.49198, birthDeathModelSubCrit.calculateLogNodeAgeDensity(5), 1e-4,
                "Log density incorrect.");
        assertEquals(-3.425174, birthDeathModelSubCrit.calculateLogNodeAgeDensity(1), 1e-4,
                "Log density incorrect.");

        assertEquals(-4.795791, birthDeathModelCrit.calculateLogNodeAgeDensity(10), 1e-4,
                "Log density incorrect.");
        assertEquals(-3.583519, birthDeathModelCrit.calculateLogNodeAgeDensity(5), 1e-4,
                "Log density incorrect.");
        assertEquals(-1.386294, birthDeathModelCrit.calculateLogNodeAgeDensity(1), 1e-4,
                "Log density incorrect.");
    }
    
    @Test
    public void calculateLogNodeAgeDensityTest() {
        assertEquals(0.0, birthDeathModelSuperCrit.calculateLogNodeAgeCDF(1e8), 1e-6,
                "Log CDF incorrect.");
        assertEquals(Double.NEGATIVE_INFINITY, birthDeathModelSuperCrit.calculateLogNodeAgeCDF(0.0), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.0002269842, birthDeathModelSuperCrit.calculateLogNodeAgeCDF(10), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.03335573, birthDeathModelSuperCrit.calculateLogNodeAgeCDF(5), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-1.363508, birthDeathModelSuperCrit.calculateLogNodeAgeCDF(1), 1e-6,
                "Log CDF incorrect.");

        assertEquals(-2.3978952727983707, birthDeathModelSubCrit.calculateLogNodeAgeCDF(1e8), 1e-6,
                "Log CDF incorrect.");
        assertEquals(Double.NEGATIVE_INFINITY, birthDeathModelSubCrit.calculateLogNodeAgeCDF(0.0), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-2.397937, birthDeathModelSubCrit.calculateLogNodeAgeCDF(10), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-2.404043, birthDeathModelSubCrit.calculateLogNodeAgeCDF(5), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-2.822555, birthDeathModelSubCrit.calculateLogNodeAgeCDF(1), 1e-6,
                "Log CDF incorrect.");

        assertEquals(0.0, birthDeathModelCrit.calculateLogNodeAgeCDF(1e8), 1e-6,
                "Log CDF incorrect.");
        assertEquals(Double.NEGATIVE_INFINITY, birthDeathModelCrit.calculateLogNodeAgeCDF(0.0), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.09531018, birthDeathModelCrit.calculateLogNodeAgeCDF(10), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.1823216, birthDeathModelCrit.calculateLogNodeAgeCDF(5), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.6931472, birthDeathModelCrit.calculateLogNodeAgeCDF(1), 1e-6,
                "Log CDF incorrect.");
    }
}
