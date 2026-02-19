package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

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
   public void calculateTreeLogLikelihoodTest() {
      RealParameter birthRate = new RealParameter("2.0");
      RealParameter deathRate = new RealParameter("1.0");
      RealParameter rhoVal = new RealParameter("0.5");

      SkylineParameter birthRateSkyline = createSkyline(birthRate);
      SkylineParameter deathRateSkyline = createSkyline(deathRate);

      // --- PART 1: Success Case ---
      CalibratedBirthDeathSkylineModel bdskyRoot = new CalibratedBirthDeathSkylineModel();
      bdskyRoot.initByName(
              "birthRate", birthRateSkyline,
              "deathRate", deathRateSkyline,
              "rho", rhoVal,
              "tree", tree,
              "conditionOnRoot", true
      );

      CalibratedBirthDeathModel constantBD = new CalibratedBirthDeathModel();
      constantBD.initByName("birthRate", birthRate, "deathRate", deathRate,
              "rho", rhoVal, "tree", tree, "conditionOnRoot", true);

      assertEquals(constantBD.calculateTreeLogLikelihood(tree),
              bdskyRoot.calculateTreeLogLikelihood(tree), 1e-50);

      // --- PART 2: Failure Case (Invalid Origin) ---
      CalibratedBirthDeathSkylineModel bdskyInvalidOrigin = new CalibratedBirthDeathSkylineModel();
      bdskyInvalidOrigin.initByName(
              "birthRate", birthRateSkyline,
              "deathRate", deathRateSkyline,
              "rho", rhoVal,
              "tree", tree,
              "origin", new RealParameter("3.5"), // < Tree Height
              "conditionOnRoot", false
      );

      // This returns -Infinity immediately
      bdskyInvalidOrigin.originInput.setValue(new RealParameter("1.5"), bdskyInvalidOrigin);
      double val = bdskyInvalidOrigin.calculateTreeLogLikelihood(tree);
      assertEquals(Double.NEGATIVE_INFINITY, val, "Should return -Infinity if Origin < Tree Height");


      // --- PART 3: Dynamic Updates ---
      RealParameter safeOrigin = new RealParameter("5.0");
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
      model.initByName(
              "birthRate", createSkyline(b),
              "deathRate", createSkyline(d),
              "rho", rho,
              "tree", tree,
              "origin", safeOrigin
      );

      Tree tallTree = new TreeParser();
      tallTree.initByName("newick", "((A:3,B:3):3,C:6);",
              "IsLabelledNewick", true, "adjustTipHeights", true);

      // Update tree
      model.treeInput.setValue(tallTree, model);

      double logL = model.calculateTreeLogLikelihood(tallTree);
      assertEquals(Double.NEGATIVE_INFINITY, logL);

      // Fix origin
      model.originInput.setValue(new RealParameter("6.1"), model);

      // Update rates
      RealParameter newRates = new RealParameter("3.0 2.0");
      RealParameter newTimes = new RealParameter("0.5");
      SkylineParameter newBirthSkyline = createSkyline(newRates, newTimes, false, true);

      model.birthRateInput.setValue(newBirthSkyline, model);

      double result = model.calculateTreeLogLikelihood(tallTree);
      assertTrue(Double.isFinite(result));

      CalibratedBirthDeathSkylineModel bdskyModel = new CalibratedBirthDeathSkylineModel();
      bdskyModel.initByName("birthRate", createSkyline(new RealParameter("3.0 2.0")),
              "deathRate", createSkyline(new RealParameter("1.0")),
              "rho", new RealParameter("0.5"),
              "tree", tallTree,
              "origin", new RealParameter("6.1"));
      System.out.println("logP = " + bdskyModel.calculateTreeLogLikelihood(tallTree));
   }
}