package calibratedcpp;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.speciation.CalibratedBirthDeathModel;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.model.BirthDeathModel;
import calibration.CalibrationClade;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

@State(Scope.Thread) // Each thread gets its own instance of the state
@BenchmarkMode(Mode.AverageTime) // Measure average time per operation
@OutputTimeUnit(TimeUnit.MICROSECONDS) // Output results in microseconds (us)
@Warmup(iterations = 3, time = 1) // 3 warm-up iterations of 1 second each
@Measurement(iterations = 5, time = 1) // 5 measurement iterations of 1 second each
@Fork(1) // Run in a separate JVM process to ensure clean environment
public class LikelihoodBenchmark {

    @Param({"10"})
    public int nCalibrations;

    CalibratedCoalescentPointProcess cpp;
    BirthDeathModel birthDeathModel;
    CalibratedBirthDeathModel heled_and_drummond;
    Tree tree;
    List<CalibrationClade> calibrationsClades;
    List<CalibrationPoint> calibrationPoints;

    // This runs once before the benchmarks start
    // We do initialization here so it doesn't count towards the execution time
    @Setup(Level.Trial)
    public void setup() {
        cpp = new CalibratedCoalescentPointProcess();
        heled_and_drummond = new CalibratedBirthDeathModel();
        birthDeathModel = new BirthDeathModel();
        calibrationsClades = new ArrayList<>();
        calibrationPoints = new ArrayList<>();
        tree = new TreeParser();
        String newick = "(((((((leaf_1:1,leaf_2:1):1,(leaf_3:1,leaf_4:1):1):1,((leaf_5:1,leaf_6:1):1,(leaf_7:1,leaf_8:1):1):1):1,(((leaf_9:1,leaf_10:1):1,(leaf_11:1,leaf_12:1):1):1,((leaf_13:1,leaf_14:1):1,(leaf_15:1,leaf_16:1):1):1):1):1,((((leaf_17:1,leaf_18:1):1,(leaf_19:1,leaf_20:1):1):1,((leaf_21:1,leaf_22:1):1,(leaf_23:1,leaf_24:1):1):1):1,(((leaf_25:1,leaf_26:1):1,(leaf_27:1,leaf_28:1):1):1,((leaf_29:1,leaf_30:1):1,(leaf_31:1,leaf_32:1):1):1):1):1):1,(((((leaf_33:1,leaf_34:1):1,(leaf_35:1,leaf_36:1):1):1,((leaf_37:1,leaf_38:1):1,(leaf_39:1,leaf_40:1):1):1):1,(((leaf_41:1,leaf_42:1):1,(leaf_43:1,leaf_44:1):1):1,((leaf_45:1,leaf_46:1):1,(leaf_47:1,leaf_48:1):1):1):1):1,((((leaf_49:1,leaf_50:1):1,(leaf_51:1,leaf_52:1):1):1,((leaf_53:1,leaf_54:1):1,(leaf_55:1,leaf_56:1):1):1):1,(((leaf_57:1,leaf_58:1):1,(leaf_59:1,leaf_60:1):1):1,((leaf_61:1,leaf_62:1):1,(leaf_63:1,leaf_64:1):1):1):1):1):1):1,((((((leaf_65:1,leaf_66:1):1,(leaf_67:1,leaf_68:1):1):1,((leaf_69:1,leaf_70:1):1,(leaf_71:1,leaf_72:1):1):1):1,(((leaf_73:1,leaf_74:1):1,(leaf_75:1,leaf_76:1):1):1,((leaf_77:1,leaf_78:1):1,(leaf_79:1,leaf_80:1):1):1):1):1,((((leaf_81:1,leaf_82:1):1,(leaf_83:1,leaf_84:1):1):1,((leaf_85:1,leaf_86:1):1,(leaf_87:1,leaf_88:1):1):1):1,(((leaf_89:1,leaf_90:1):1,(leaf_91:1,leaf_92:1):1):1,((leaf_93:1,leaf_94:1):1,(leaf_95:1,leaf_96:1):1):1):1):1):1,(((((leaf_97:1,leaf_98:1):1,(leaf_99:1,leaf_100:1):1):1,((leaf_101:1,leaf_102:1):1,(leaf_103:1,leaf_104:1):1):1):1,(((leaf_105:1,leaf_106:1):1,(leaf_107:1,leaf_108:1):1):1,((leaf_109:1,leaf_110:1):1,(leaf_111:1,leaf_112:1):1):1):1):1,((((leaf_113:1,leaf_114:1):1,(leaf_115:1,leaf_116:1):1):1,((leaf_117:1,leaf_118:1):1,(leaf_119:1,leaf_120:1):1):1):1,(((leaf_121:1,leaf_122:1):1,(leaf_123:1,leaf_124:1):1):1,((leaf_125:1,leaf_126:1):1,(leaf_127:1,leaf_128:1):1):1):1):1):1):1);";
        tree.initByName("newick", newick,
                "IsLabelledNewick", true,
                "adjustTipHeights", false);

        TaxonSet taxa = tree.getTaxonset();
        ParametricDistribution ageDist = new Uniform();
        ageDist.initByName("upper", 1.5, "lower", 0.5);

        CalibrationPoint rootCalibrationPoint = new CalibrationPoint();
        ParametricDistribution rootAgeDist = new Uniform();
        rootAgeDist.initByName("lower", 6.5, "upper", 7.5);
        rootCalibrationPoint.initByName("taxonset", taxa,
                "distr", rootAgeDist);
        calibrationPoints.add(rootCalibrationPoint);

        CalibrationPoint specialCalibrationPoint = new CalibrationPoint();
        CalibrationClade specialCalibrationClade = new CalibrationClade();

        TaxonSet specialTaxonSet = new TaxonSet();
        List<Taxon> specialTaxonList = new ArrayList<>();
        for (int i = 1; i <= 16; i++) {
            String taxon = "leaf_" + i;
            specialTaxonList.add(taxa.getTaxon(taxon));
        }
        specialTaxonSet.initByName("taxon", specialTaxonList);
        ParametricDistribution specialDist = new Uniform();
        specialDist.initByName("upper", 4.5, "lower", 3.5);
        specialCalibrationPoint.initByName("taxonset", specialTaxonSet,
                "distr", specialDist);
        specialCalibrationClade.initByName("taxa", specialTaxonSet);

        calibrationPoints.add(specialCalibrationPoint);
        calibrationsClades.add(specialCalibrationClade);

        for (int i = 0; i < nCalibrations ; i++) {
            String t1Name = "leaf_" + (2 * i + 1);
            String t2Name = "leaf_" + (2 * (i + 1));

            TaxonSet taxonSet = new TaxonSet();
            List<Taxon> taxonList = new ArrayList<>();

            taxonList.add(taxa.getTaxon(t1Name));
            taxonList.add(taxa.getTaxon(t2Name));

            taxonSet.initByName("taxon", taxonList);

            CalibrationClade calibrationClade = new CalibrationClade();
            calibrationClade.initByName("taxa", taxonSet);
            calibrationsClades.add(calibrationClade);

            CalibrationPoint calibrationPoint = new CalibrationPoint();
            calibrationPoint.initByName("taxonset", taxonSet, "distr", ageDist);
            calibrationPoints.add(calibrationPoint);
        }

        RealParameter turnover = new RealParameter("1.1");
        RealParameter birthRate = new RealParameter("2.0");
        RealParameter rho = new RealParameter("1.0");

        birthDeathModel.initByName("birthRate", birthRate,
                "turnover", turnover,
                "rho", rho);

        cpp.initByName("tree", tree,
                "treeModel", birthDeathModel,
                "calibrations", calibrationsClades,
                "conditionOnRoot", true);
        heled_and_drummond.initByName("tree", tree,
                "calibrations", calibrationPoints,
                "birthRate", birthRate,
                "relativeDeathRate", turnover,
                "sampleProbability", rho);

        double cppVal = cpp.calculateTreeLogLikelihood(tree);
        double hdVal = heled_and_drummond.calculateTreeLogLikelihood(tree);
        double epsilon = 1e-9;

        double logNfactorial = 0.0;
        for (int i = 1; i <= taxa.getTaxonCount() ; i++) {
            logNfactorial += Math.log(i);
        }

        if (Math.abs(cppVal - logNfactorial - hdVal) > epsilon) {
            throw new RuntimeException(String.format(
                    "Likelihood mismatch! CPP: %f, BD: %f", cppVal, hdVal
            ));
        }

        System.out.println("Verification passed: Both models produce " + cppVal);
    }

    // Benchmark for CalibratedBirthDeathModel
    @Benchmark
    public double measureHeledAndDrummond() {
        return heled_and_drummond.calculateTreeLogLikelihood(tree);
    }

    // Benchmark for CalibratedCoalescentPointProcess
    @Benchmark
    public double measureCPP() {
        return cpp.calculateTreeLogLikelihood(tree);
    }

    // Standard main method to launch the benchmark
    public static void main(String[] args) throws Exception {
        // -------------------------------------------------------------
        // 3. GENERATE PARAMETERS PROGRAMMATICALLY
        // -------------------------------------------------------------
        // Create an array of strings: ["1", "2", "3", ... "50"]
        String[] calibrationCounts = IntStream.rangeClosed(1, 2)
                .mapToObj(String::valueOf)
                .toArray(String[]::new);

        Options opt = new OptionsBuilder()
                .include(LikelihoodBenchmark.class.getSimpleName())
                // Pass the sequence 1..50 to the "nCalibrations" field
                .param("nCalibrations", calibrationCounts)
                .build();

        new Runner(opt).run();
    }
}