package calibration;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import calibrationprior.CalibrationCladePrior;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.*;

public class CalibrationForestParserTest {

    private static CalibrationForestParser parse(String newick) {
        CalibrationForestParser ct = new CalibrationForestParser();
        ct.initByName("newick", newick);
        return ct;
    }

    private static Set<String> names(TaxonSet ts) {
        return ts.getTaxonSet().stream().map(Taxon::getID).collect(Collectors.toSet());
    }

    // ==================== TaxonSet output ====================

    @Test
    void singleCalibration_taxonSet() {
        CalibrationForestParser ct = parse("(A,B)[&name=Root,lower=10.0,upper=20.0];");

        List<TaxonSet> sets = ct.getTaxonSets();
        assertEquals(1, sets.size());
        assertEquals("Root", sets.get(0).getID());
        assertEquals(Set.of("A", "B"), names(sets.get(0)));
        // hasVirtualRoot() was removed — the root in this case is the single calibration node, not a virtual root
    }

    @Test
    void twoNestedCalibrations_taxonSets() {
        CalibrationForestParser ct = parse("((A,B)[&name=Inner,lower=5.0,upper=10.0],C)[&name=Outer,lower=15.0,upper=30.0];");

        List<TaxonSet> sets = ct.getTaxonSets();
        assertEquals(2, sets.size());

        TaxonSet inner = sets.stream().filter(s -> "Inner".equals(s.getID())).findFirst().orElseThrow();
        TaxonSet outer = sets.stream().filter(s -> "Outer".equals(s.getID())).findFirst().orElseThrow();

        assertEquals(Set.of("A", "B"), names(inner));
        assertEquals(Set.of("A", "B", "C"), names(outer));
    }

    // ==================== CalibrationCladePrior output ====================

    @Test
    void boundsProduceCalibrationCladePrior() {
        CalibrationForestParser ct = parse("(A,B)[&name=Root,lower=10.0,upper=20.0];");

        List<CalibrationCladePrior> priors = ct.getCalibrationCladePriors();
        assertEquals(1, priors.size());
        assertEquals("Root", priors.get(0).getID());
        assertEquals(10.0, priors.get(0).getLower(), 1e-9);
        assertEquals(20.0, priors.get(0).getUpper(), 1e-9);
    }

    @Test
    void nodeWithoutBounds_noCalibrationCladePrior() {
        // Root has neither name nor bounds — no TaxonSet, no prior
        CalibrationForestParser ct = parse("((A,B)[&name=C1,lower=1.0,upper=2.0],(C,D)[&name=C2,lower=3.0,upper=4.0]);");

        assertEquals(2, ct.getTaxonSets().size());
        assertEquals(2, ct.getCalibrationCladePriors().size());
    }

    @Test
    void namedNodeWithoutBounds_taxonSetOnly() {
        // Has name but no bounds — gets a TaxonSet but no CalibrationCladePrior
        CalibrationForestParser ct = parse("(A,B)[&name=TopClade];");

        assertEquals(1, ct.getTaxonSets().size());
        assertEquals(0, ct.getCalibrationCladePriors().size());
        assertEquals("TopClade", ct.getTaxonSets().get(0).getID());
    }

    // ==================== Virtual root ====================

    @Test
    void virtualRoot_excludedFromBothOutputs() {
        CalibrationForestParser ct = parse(
                "((A,B)[&name=Clade1,lower=10.0,upper=20.0],(C,D)[&name=Clade2,lower=5.0,upper=15.0])[&virtualRoot=true];");

        // Virtual root is excluded from taxon sets and priors
        assertEquals(2, ct.getTaxonSets().size());
        assertEquals(2, ct.getCalibrationCladePriors().size());
        assertEquals(Set.of("A", "B", "C", "D"), ct.getAllTaxa());
    }

    // ==================== Alternative key spelling ====================

    @Test
    void lowerUpperAgeVariants() {
        CalibrationForestParser ct = parse("(A,B)[&name=X,lowerAge=8.0,upperAge=16.0];");

        CalibrationCladePrior prior = ct.getCalibrationCladePriors().get(0);
        assertEquals(8.0, prior.getLower(), 1e-9);
        assertEquals(16.0, prior.getUpper(), 1e-9);
    }

    // ==================== Comment stripping ====================

    @Test
    void commentLinesIgnored() {
        String newick = "# comment\n(A,B)[&name=X,lower=1.0,upper=2.0];";
        CalibrationForestParser ct = parse(newick);
        assertEquals(1, ct.getTaxonSets().size());
    }

    // ==================== Branch lengths ====================

    @Test
    void branchLengthsIgnored() {
        CalibrationForestParser ct = parse("((A:1.0,B:2.0):0.5)[&name=Cl,lower=5.0,upper=10.0];");
        assertEquals(Set.of("A", "B"), names(ct.getTaxonSets().get(0)));
    }

    // ==================== All taxa ====================

    @Test
    void allTaxaCollected() {
        CalibrationForestParser ct = parse("((A,B)[&lower=1.0,upper=2.0],(C,D,E));");
        assertEquals(Set.of("A", "B", "C", "D", "E"), ct.getAllTaxa());
    }

    // ==================== TaxonSet and prior share the same TaxonSet instance ====================

    @Test
    void priorSharesTaxonSetWithTaxonSetList() {
        CalibrationForestParser ct = parse("(A,B)[&name=Root,lower=1.0,upper=2.0];");

        TaxonSet fromSets = ct.getTaxonSets().get(0);
        TaxonSet fromPrior = ct.getCalibrationCladePriors().get(0).getTaxa();

        assertSame(fromSets, fromPrior);
    }
}
