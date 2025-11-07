package calibratedcpp.model;

import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class BirthDeathModelTest {

    private final BirthDeathModel birthDeathModelSuperCrit;
    private final BirthDeathModel birthDeathModelSubCrit;
    private final BirthDeathModel birthDeathModelCrit;

    public BirthDeathModelTest() {
        birthDeathModelSuperCrit = new BirthDeathModel();
        birthDeathModelSuperCrit.initByName("birthRate", new RealParameter("2.0"), "deathRate", new RealParameter("1.0"), "rho", new RealParameter("0.1"));

        birthDeathModelSubCrit = new BirthDeathModel();
        birthDeathModelSubCrit.initByName("birthRate", new RealParameter("1.0"), "deathRate", new RealParameter("2.0"), "rho", new RealParameter("0.1"));

        birthDeathModelCrit = new BirthDeathModel();
        birthDeathModelCrit.initByName("birthRate", new RealParameter("1.0"), "diversificationRate", new RealParameter("0.0"), "rho", new RealParameter("1.0"));
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
    public void calculateLogDensity() {
        assertEquals(-8.390925, birthDeathModelSuperCrit.calculateLogDensity(10), 1e-4,
                "Log density incorrect.");
        assertEquals(-3.443752, birthDeathModelSuperCrit.calculateLogDensity(5), 1e-4,
                "Log density incorrect.");
        assertEquals(-1.200227, birthDeathModelSuperCrit.calculateLogDensity(1), 1e-4,
                "Log density incorrect.");

        assertEquals(-12.4932, birthDeathModelSubCrit.calculateLogDensity(10), 1e-4,
                "Log density incorrect.");
        assertEquals(-7.49198, birthDeathModelSubCrit.calculateLogDensity(5), 1e-4,
                "Log density incorrect.");
        assertEquals(-3.425174, birthDeathModelSubCrit.calculateLogDensity(1), 1e-4,
                "Log density incorrect.");

        assertEquals(-4.795791, birthDeathModelCrit.calculateLogDensity(10), 1e-4,
                "Log density incorrect.");
        assertEquals(-3.583519, birthDeathModelCrit.calculateLogDensity(5), 1e-4,
                "Log density incorrect.");
        assertEquals(-1.386294, birthDeathModelCrit.calculateLogDensity(1), 1e-4,
                "Log density incorrect.");
    }

    @Test
    public void calculateLogCDF() {
        assertEquals(0.0, birthDeathModelSuperCrit.calculateLogCDF(1e8), 1e-6,
                "Log CDF incorrect.");
        assertEquals(Double.NEGATIVE_INFINITY, birthDeathModelSuperCrit.calculateLogCDF(0.0), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.0002269842, birthDeathModelSuperCrit.calculateLogCDF(10), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.03335573, birthDeathModelSuperCrit.calculateLogCDF(5), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-1.363508, birthDeathModelSuperCrit.calculateLogCDF(1), 1e-6,
                "Log CDF incorrect.");

        assertEquals(-2.3978952727983707, birthDeathModelSubCrit.calculateLogCDF(1e8), 1e-6,
                "Log CDF incorrect.");
        assertEquals(Double.NEGATIVE_INFINITY, birthDeathModelSubCrit.calculateLogCDF(0.0), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-2.397937, birthDeathModelSubCrit.calculateLogCDF(10), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-2.404043, birthDeathModelSubCrit.calculateLogCDF(5), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-2.822555, birthDeathModelSubCrit.calculateLogCDF(1), 1e-6,
                "Log CDF incorrect.");

        assertEquals(0.0, birthDeathModelCrit.calculateLogCDF(1e8), 1e-6,
                "Log CDF incorrect.");
        assertEquals(Double.NEGATIVE_INFINITY, birthDeathModelCrit.calculateLogCDF(0.0), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.09531018, birthDeathModelCrit.calculateLogCDF(10), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.1823216, birthDeathModelCrit.calculateLogCDF(5), 1e-6,
                "Log CDF incorrect.");
        assertEquals(-0.6931472, birthDeathModelCrit.calculateLogCDF(1), 1e-6,
                "Log CDF incorrect.");
    }
}