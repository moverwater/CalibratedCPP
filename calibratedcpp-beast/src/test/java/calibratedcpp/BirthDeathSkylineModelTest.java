package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class BirthDeathSkylineModelTest {

   private Tree tree;
   private RealParameter b, d, rho, origin;

   @BeforeEach
   public void setUp() {
      // 1. Initialize a simple fixed tree
      tree = new TreeParser();
      tree.initByName("newick", "((A:1,B:1):1,C:2);",
              "IsLabelledNewick", true,
              "adjustTipHeights", true);

      // 2. Initialize standard parameters
      b = new RealParameter("2.0");
      d = new RealParameter("1.0");
      rho = new RealParameter("1.0");
      origin = new RealParameter("5.0");
   }

   @Test
   public void testOnlyOneRateThrows() {
      BirthDeathSkylineModel model = new BirthDeathSkylineModel();

      // We only provide ONE rate (Birth), missing the second required rate
      model.setInputValue("birthRate", b);
      model.setInputValue("rho", rho);
      model.setInputValue("tree", tree);
      model.setInputValue("origin", origin); // Always provide origin to satisfy parent class
      model.setInputValue("conditionOnRoot", true);

      // Expect validation failure
      assertThrows(IllegalArgumentException.class, model::initAndValidate,
              "Should throw because exactly TWO rates must be specified.");
   }

   @Test
   public void testThreeRatesThrows() {
      BirthDeathSkylineModel model = new BirthDeathSkylineModel();
      RealParameter div = new RealParameter("0.5");

      // We provide THREE rates (Birth, Death, Diversification)
      model.setInputValue("birthRate", b);
      model.setInputValue("deathRate", d);
      model.setInputValue("diversificationRate", div);
      model.setInputValue("rho", rho);
      model.setInputValue("tree", tree);
      model.setInputValue("origin", origin);
      model.setInputValue("conditionOnRoot", true);

      assertThrows(IllegalArgumentException.class, model::initAndValidate,
              "Should throw because mutual exclusion check failed (3 rates provided).");
   }

   @Test
   public void testRelativeOutOfBoundsThrows() {
      BirthDeathSkylineModel model = new BirthDeathSkylineModel();

      // 1.5 is invalid for a relative time (must be 0.0 to 1.0)
      RealParameter badTime = new RealParameter("1.5");

      model.setInputValue("birthRate", b);
      model.setInputValue("deathRate", d);

      // Set the bad time and flag it as relative
      model.setInputValue("birthRateChangeTimes", badTime);
      model.setInputValue("birthRateTimesRelative", true);

      model.setInputValue("rho", rho);
      model.setInputValue("tree", tree);
      model.setInputValue("origin", origin);
      model.setInputValue("conditionOnRoot", true);

      assertThrows(IllegalArgumentException.class, model::initAndValidate,
              "Should throw because relative time 1.5 is outside [0, 1].");
   }

   @Test
   public void testLikelihoodMatchesConstantBD() {
      // 1. Setup BDSKY (this model)
      BirthDeathSkylineModel bdsky = new BirthDeathSkylineModel();
      bdsky.initByName("birthRate", b, "deathRate", d, "rho", rho,
              "tree", tree, "origin", origin, "conditionOnRoot", true);

      // 2. Setup Reference Model (Standard Birth-Death)
      // Ensure BirthDeathModel is available in your classpath
      BirthDeathModel constantBD = new BirthDeathModel();
      constantBD.initByName("birthRate", b, "deathRate", d, "rho", rho,
              "tree", tree, "origin", origin, "conditionOnRoot", true);

      // 3. Compare Likelihoods
      // Because BDSKY with 0 change times IS a constant Birth-Death model,
      // the likelihoods must match to high precision.
      double logL_BDSKY = bdsky.calculateTreeLogLikelihood(tree);
      double logL_BD = constantBD.calculateTreeLogLikelihood(tree);

      assertEquals(logL_BD, logL_BDSKY, 1e-10,
              "BDSKY likelihood should match constant Birth-Death when no intervals are defined.");
   }

   @Test
   public void calculateLogNodeAgeDensityTest() {
      BirthDeathSkylineModel model = new BirthDeathSkylineModel();
      model.initByName("birthRate", b, "deathRate", d, "rho", rho,
              "tree", tree, "origin", origin, "conditionOnRoot", true);

      // Force a calculation to populate internal arrays
      model.calculateTreeLogLikelihood(tree);

      // Check density at an arbitrary time t=0.5
      double density = model.calculateLogNodeAgeDensity(0.5);

      assertTrue(Double.isFinite(density), "Density should be a finite log value.");
   }

   @Test
   public void calculateLogNodeAgeCDFTest() {
      BirthDeathSkylineModel model = new BirthDeathSkylineModel();
      model.initByName("birthRate", b, "deathRate", d, "rho", rho,
              "tree", tree, "origin", origin, "conditionOnRoot", true);

      model.calculateTreeLogLikelihood(tree);

      // At present (t=0), no lineages have formed yet in the backward view, CDF = 0, log(CDF) = -Inf
      double cdfZero = model.calculateLogNodeAgeCDF(0.0);
      assertEquals(Double.NEGATIVE_INFINITY, cdfZero, 1e-10, "log(CDF) at t=0 should be -Infinity.");

      // At the root, CDF should be <= 1.0 (log <= 0)
      double cdfRoot = model.calculateLogNodeAgeCDF(tree.getRoot().getHeight());
      assertTrue(cdfRoot <= 0.0, "log(CDF) should never be positive.");
   }

   @Test
   public void calculateTreeLogLikelihoodTest() {
      RealParameter birthRate = new RealParameter("2.0");
      RealParameter deathRate = new RealParameter("1.0");
      RealParameter rho = new RealParameter("0.5");

      // --- PART 1: Success Case (Comparison) ---
      BirthDeathSkylineModel bdskyRoot = new BirthDeathSkylineModel();
      bdskyRoot.initByName("birthRate", birthRate, "deathRate", deathRate,
              "rho", rho, "tree", tree, "conditionOnRoot", true);

      BirthDeathModel constantBD = new BirthDeathModel();
      constantBD.initByName("birthRate", birthRate, "deathRate", deathRate,
              "rho", rho, "tree", tree, "conditionOnRoot", true);

      assertEquals(constantBD.calculateTreeLogLikelihood(tree),
              bdskyRoot.calculateTreeLogLikelihood(tree), 1e-12,
              "BDSKY should match constant Birth-Death when no intervals are used");

      // --- PART 2: Failure Case (Invalid Origin) ---
      assertThrows(Throwable.class, () -> {
         BirthDeathSkylineModel bdskyInvalidOrigin = new BirthDeathSkylineModel();
         bdskyInvalidOrigin.initByName(
                 "birthRate", birthRate,
                 "deathRate", deathRate,
                 "rho", rho,
                 "tree", tree,
                 "origin", new RealParameter("1.5") // < Tree Height (2.0)
         );
      }, "Model should throw an exception if Origin is less than Tree Height");

      Tree validTree = new TreeParser();
      validTree.initByName("newick", "((A:1,B:1):1,C:2);",
              "IsLabelledNewick", true, "adjustTipHeights", true);

      RealParameter origin = new RealParameter("5.0"); // Origin is safely larger

      BirthDeathSkylineModel model = new BirthDeathSkylineModel();
      model.initByName(
              "birthRate", b,
              "deathRate", d,
              "rho", rho,
              "tree", validTree, // Pass the valid tree first
              "origin", origin,
              "conditionOnRoot", false
      );

      Tree tallTree = new TreeParser();
      tallTree.initByName("newick", "((A:3,B:3):3,C:6);",
              "IsLabelledNewick", true, "adjustTipHeights", true);

      model.treeInput.setValue(tallTree, model);

      double logL = model.calculateTreeLogLikelihood(tallTree);

      assertEquals(Double.NEGATIVE_INFINITY, logL,
              "Likelihood should be -Inf when the tree grows taller than the origin.");
   }
}