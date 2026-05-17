package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import bdsky.evolution.speciation.BirthDeathSkylineModel;
import static org.junit.jupiter.api.Assertions.*;

public class CalibratedBirthDeathSkylineModelTest {

   private Tree tree;
   private RealParameter b, d, rho, origin;

   @BeforeEach
   public void setUp() {
      tree = new TreeParser();
      tree.initByName("newick", "((A:1,B:1):1,C:2);",
              "IsLabelledNewick", true,
              "adjustTipHeights", true);

      b = new RealParameter("2.0");
      d = new RealParameter("1.0");
      rho = new RealParameter("1.0");
      origin = new RealParameter("5.0");
   }

   private SkylineParameter createSkyline(RealParameter rates, RealParameter times, boolean relative, boolean reverse) {
      SkylineParameter sp = new SkylineParameter();
      sp.initByName("values", rates,
              "changeTimes", times,
              "timesAreRelative", relative,
              "timesAreAges", reverse);
      return sp;
   }

   private SkylineParameter createSkyline(RealParameter rates) {
      return createSkyline(rates, null, false, false);
   }

   @Test
   public void testOnlyOneRateThrows() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();

       assertThrows(RuntimeException.class, () -> model.initByName(
               "birthRate", createSkyline(b),
               "rho", rho,
               "tree", tree,
               "origin", origin,
               "conditionOnRoot", true
       ), "Should throw because exactly TWO rates must be specified.");

      // assertEquals(IllegalArgumentException.class, exception.getCause().getClass());
   }

   @Test
   public void testThreeRatesThrows() {
      RealParameter div = new RealParameter("0.5");
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();

      assertThrows(RuntimeException.class, () -> model.initByName(
              "birthRate", createSkyline(b),
              "deathRate", createSkyline(d),
              "diversificationRate", createSkyline(div),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      ), "Should throw because mutual exclusion check failed (3 rates provided).");
   }

   @Test
   public void testRelativeOutOfBoundsThrows() {
      RealParameter badTime = new RealParameter("1.5");
      RealParameter rates = new RealParameter("2.0 1.0");

      assertThrows(RuntimeException.class, () -> createSkyline(rates, badTime, true, false), "Should throw because relative time 1.5 is outside [0, 1].");
   }

   @Test
   public void testLikelihoodMatchesConstantBD() {
      CalibratedBirthDeathSkylineModel bdsky = new CalibratedBirthDeathSkylineModel();
      bdsky.initByName(
              "birthRate", createSkyline(b),
              "deathRate", createSkyline(d),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      );

      CalibratedBirthDeathModel constantBD = new CalibratedBirthDeathModel();
      constantBD.initByName("birthRate", b, "deathRate", d, "rho", rho,
              "tree", tree, "origin", origin, "conditionOnRoot", true);

      double logL_BDSKY = bdsky.calculateTreeLogLikelihood(tree);
      double logL_BD = constantBD.calculateTreeLogLikelihood(tree);

      assertEquals(logL_BD, logL_BDSKY, 1e-10,
              "BDSKY likelihood should match constant Birth-Death when no intervals are defined.");
   }

   @Test
   public void calculateLogNodeAgeDensityTest() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
      model.initByName(
              "birthRate", createSkyline(b),
              "deathRate", createSkyline(d),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      );

      model.calculateTreeLogLikelihood(tree);
      double density = model.calculateLogNodeAgeDensity(0.5);
      assertTrue(Double.isFinite(density), "Density should be a finite log value.");
   }

   @Test
   public void calculateLogNodeAgeCDFTest() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
      model.initByName(
              "birthRate", createSkyline(b),
              "deathRate", createSkyline(d),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      );

      model.calculateTreeLogLikelihood(tree);

      double cdfZero = model.calculateLogNodeAgeCDF(0.0);
      assertEquals(Double.NEGATIVE_INFINITY, cdfZero, 1e-10, "log(CDF) at t=0 should be -Infinity.");

      double cdfRoot = model.calculateLogNodeAgeCDF(tree.getRoot().getHeight());
      assertTrue(cdfRoot <= 0.0, "log(CDF) should never be positive.");
   }

   @Test
   public void calculateTreeLogLikelihood() {
      RealParameter birthRates = new RealParameter("2.0 1.0 3.0");
      RealParameter deathRates = new RealParameter("1.1 2.0 0.5");
      RealParameter rho = new RealParameter("0.5");
      TreeInterface tree = new TreeParser("((A:3,B:3):4,(C:6,D:6):1);");

      RealParameter birthRateChangeTimes = new RealParameter("1.0 1.5");
      RealParameter deathRateChangeTimes = new RealParameter("0.5 1.25");

      RealParameter birthRateChangeTimesBDSKY = new RealParameter("0.0 1.0 1.5");
      RealParameter deathRateChangeTimesBDSKY = new RealParameter("0.0 0.5 1.25");

      BirthDeathSkylineModel BDSKY = new BirthDeathSkylineModel();
      BDSKY.initByName("birthRate", birthRates,
              "deathRate", deathRates,
              "rho", rho,
              "birthRateChangeTimes", birthRateChangeTimesBDSKY,
              "deathRateChangeTimes", deathRateChangeTimesBDSKY,
              "samplingRate", new RealParameter("0.0"),
              "origin", new RealParameter("8.0"),
              "conditionOnSurvival", true,
              "reverseTimeArrays", new BooleanParameter("true true true true true"),
              "tree", tree);

      CalibratedBirthDeathSkylineModel CalibratedBDSKY = new CalibratedBirthDeathSkylineModel();
      SkylineParameter birthRateSkyline = new SkylineParameter();
      SkylineParameter deathRateSkyline = new SkylineParameter();
      birthRateSkyline.initByName("values", birthRates,
              "changeTimes", birthRateChangeTimes,
              "timesAreAges", true);
      deathRateSkyline.initByName("values", deathRates,
              "changeTimes", deathRateChangeTimes,
              "timesAreAges", true);
      CalibratedBDSKY.initByName("birthRate", birthRateSkyline,
              "deathRate", deathRateSkyline,
              "rho", rho,
              "tree", tree,
              "origin", new RealParameter("8.0"));

      assertEquals(BDSKY.calculateTreeLogLikelihood(tree) + 3.0 * Math.log(2.0) - Math.log(4.0) - Math.log(3.0) - Math.log(2.0),
              CalibratedBDSKY.calculateTreeLogLikelihood(tree), 1e-8,
              "Likelihood of tree under BDSKY does not match likelihood under Calibrated BDSKY.");
   }
}