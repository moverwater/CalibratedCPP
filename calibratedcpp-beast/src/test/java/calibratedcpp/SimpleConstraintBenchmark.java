package calibratedcpp;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import calibratedcpp.model.BirthDeathModel;
import calibration.CalibrationClade;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.regex.*;

/**
 * Simple benchmark that reads calibration constraints from annotated Newick files.
 * No JMH dependency - just basic timing.
 *
 * Default: reads from validation/phylodata_calibration_forests/
 *
 * Usage:
 *   java -cp "lib/*:target/classes" calibratedcpp.SimpleConstraintBenchmark [newick-directory] [iterations] [warmup]
 */
public class SimpleConstraintBenchmark {

    private static final String DEFAULT_CONSTRAINTS_DIR = "validation/phylodata_calibration_forests";
    private static final Pattern ANNOTATION_PATTERN = Pattern.compile("\\[&([^\\]]+)\\]");
    private static final Pattern KEY_VALUE_PATTERN = Pattern.compile("([^=,]+)=([^,\\]]+)");

    public static void main(String[] args) throws Exception {
        int iterations = 1000;
        int warmup = 100;
        Path inputPath = null;

        // Parse arguments - first non-numeric arg is path, rest are iterations/warmup
        int argIdx = 0;
        if (args.length > 0 && !args[0].matches("\\d+")) {
            inputPath = Path.of(args[0]);
            argIdx = 1;
        }
        if (args.length > argIdx) {
            iterations = Integer.parseInt(args[argIdx]);
        }
        if (args.length > argIdx + 1) {
            warmup = Integer.parseInt(args[argIdx + 1]);
        }

        // Use default path if not specified
        if (inputPath == null) {
            inputPath = Path.of(DEFAULT_CONSTRAINTS_DIR);
            if (!Files.exists(inputPath)) {
                inputPath = Path.of("calibratedcpp-beast", DEFAULT_CONSTRAINTS_DIR);
            }
            System.out.println("Using default constraints directory: " + inputPath);
        }
        List<Path> newickFiles = new ArrayList<>();

        if (Files.isDirectory(inputPath)) {
            try (var stream = Files.list(inputPath)) {
                stream.filter(p -> p.toString().endsWith(".newick"))
                      .sorted()
                      .forEach(newickFiles::add);
            }
        } else if (Files.exists(inputPath)) {
            newickFiles.add(inputPath);
        } else {
            System.err.println("Path not found: " + inputPath);
            return;
        }

        if (newickFiles.isEmpty()) {
            System.err.println("No .newick files found in " + inputPath);
            return;
        }

        System.out.println("Found " + newickFiles.size() + " constraint file(s)");
        System.out.println("Iterations: " + iterations + ", Warmup: " + warmup);
        System.out.println();

        // Output header
        System.out.printf("%-60s %8s %8s %12s %12s %12s%n",
                "File", "Taxa", "Calib", "Mean(us)", "Min(us)", "Max(us)");
        System.out.println("-".repeat(110));

        for (Path newickFile : newickFiles) {
            try {
                BenchmarkResult result = runBenchmark(newickFile, iterations, warmup);
                System.out.printf("%-60s %8d %8d %12.2f %12.2f %12.2f%n",
                        newickFile.getFileName(),
                        result.numTaxa,
                        result.numCalibrations,
                        result.meanTimeUs,
                        result.minTimeUs,
                        result.maxTimeUs);
            } catch (Exception e) {
                System.err.println("Error processing " + newickFile.getFileName() + ": " + e.getMessage());
            }
        }
    }

    static class BenchmarkResult {
        int numTaxa;
        int numCalibrations;
        double meanTimeUs;
        double minTimeUs;
        double maxTimeUs;
    }

    private static BenchmarkResult runBenchmark(Path newickPath, int iterations, int warmup) throws Exception {
        String content = Files.readString(newickPath);
        ParsedConstraints parsed = parseConstraintNewick(content);

        BenchmarkResult result = new BenchmarkResult();
        result.numTaxa = parsed.allTaxa.size();
        result.numCalibrations = parsed.calibrations.size();

        if (result.numTaxa == 0) {
            throw new RuntimeException("No taxa found in file");
        }

        // Create Taxon objects
        Map<String, Taxon> taxonMap = new HashMap<>();
        for (String taxonName : parsed.allTaxa) {
            Taxon taxon = new Taxon(taxonName);
            taxonMap.put(taxonName, taxon);
        }

        // Create CalibrationClades
        List<CalibrationClade> calibrationClades = new ArrayList<>();
        for (ParsedCalibration cal : parsed.calibrations) {
            if (cal.taxa.isEmpty()) continue;

            List<Taxon> taxonList = new ArrayList<>();
            for (String name : cal.taxa) {
                Taxon t = taxonMap.get(name);
                if (t != null) {
                    taxonList.add(t);
                }
            }

            if (taxonList.isEmpty()) continue;

            TaxonSet taxonSet = new TaxonSet();
            taxonSet.initByName("taxon", taxonList);

            CalibrationClade clade = new CalibrationClade();
            clade.setID(cal.name != null ? cal.name : "clade_" + calibrationClades.size());
            clade.initByName("taxa", taxonSet);
            calibrationClades.add(clade);

            // Debug: print calibration info
            // System.err.println("  Calibration: " + clade.getID() + " with " + taxonList.size() + " taxa");
        }

        if (calibrationClades.isEmpty()) {
            throw new RuntimeException("No calibrations found");
        }

        // Generate random tree
        Tree tree = generateRandomTree(parsed.allTaxa);

        // Set up birth-death model
        RealParameter turnover = new RealParameter("0.0");
        RealParameter birthRate = new RealParameter("2.0");
        RealParameter rho = new RealParameter("1.0");

        BirthDeathModel birthDeathModel = new BirthDeathModel();
        birthDeathModel.initByName("birthRate", birthRate,
                "turnover", turnover,
                "rho", rho);

        // Set up CalibratedCoalescentPointProcess
        CalibratedCoalescentPointProcess cpp = new CalibratedCoalescentPointProcess();
        cpp.initByName("tree", tree,
                "treeModel", birthDeathModel,
                "calibrations", calibrationClades,
                "conditionOnRoot", true);

        // Warmup
        for (int i = 0; i < warmup; i++) {
            cpp.calculateTreeLogLikelihood(tree);
        }

        // Benchmark
        long[] times = new long[iterations];
        for (int i = 0; i < iterations; i++) {
            long start = System.nanoTime();
            cpp.calculateTreeLogLikelihood(tree);
            times[i] = System.nanoTime() - start;
        }

        // Compute statistics
        long sum = 0;
        long min = Long.MAX_VALUE;
        long max = Long.MIN_VALUE;
        for (long t : times) {
            sum += t;
            min = Math.min(min, t);
            max = Math.max(max, t);
        }

        result.meanTimeUs = (sum / (double) iterations) / 1000.0;
        result.minTimeUs = min / 1000.0;
        result.maxTimeUs = max / 1000.0;

        return result;
    }

    private static Tree generateRandomTree(Set<String> taxaNames) throws Exception {
        List<String> taxa = new ArrayList<>(taxaNames);
        Collections.sort(taxa);

        Random rand = new Random(42);
        String newick = buildRandomNewick(taxa, rand);

        TreeParser treeParser = new TreeParser();
        treeParser.initByName("newick", newick,
                "IsLabelledNewick", true,
                "adjustTipHeights", false);

        return treeParser;
    }

    private static String buildRandomNewick(List<String> taxa, Random rand) {
        if (taxa.size() == 1) {
            return taxa.get(0) + ":0.1";
        }

        Collections.shuffle(taxa, rand);
        int splitPoint = 1 + rand.nextInt(taxa.size() - 1);

        List<String> left = new ArrayList<>(taxa.subList(0, splitPoint));
        List<String> right = new ArrayList<>(taxa.subList(splitPoint, taxa.size()));

        double branchLength = 0.1 + rand.nextDouble() * 0.5;

        String leftNewick = buildRandomNewick(left, rand);
        String rightNewick = buildRandomNewick(right, rand);

        return "(" + leftNewick + "," + rightNewick + "):" + branchLength;
    }

    // ==================== Newick Parsing ====================

    static class ParsedCalibration {
        String name;
        Set<String> taxa = new HashSet<>();
        String distType;
        Double distM, distS, distOffset;
        boolean monophyletic;
        List<ParsedCalibration> children = new ArrayList<>();
    }

    static class ParsedConstraints {
        Set<String> allTaxa = new HashSet<>();
        List<ParsedCalibration> calibrations = new ArrayList<>();
    }

    private static ParsedConstraints parseConstraintNewick(String content) {
        ParsedConstraints result = new ParsedConstraints();

        StringBuilder cleaned = new StringBuilder();
        for (String line : content.split("\n")) {
            String trimmed = line.trim();
            if (!trimmed.startsWith("#")) {
                cleaned.append(trimmed);
            }
        }

        String newick = cleaned.toString().trim();
        if (newick.endsWith(";")) {
            newick = newick.substring(0, newick.length() - 1);
        }

        if (newick.isEmpty()) {
            return result;
        }

        ParseResult parseResult = parseSubtree(newick, 0);
        if (parseResult.calibration != null) {
            collectCalibrations(parseResult.calibration, result);
        }

        return result;
    }

    private static void collectCalibrations(ParsedCalibration cal, ParsedConstraints result) {
        result.allTaxa.addAll(cal.taxa);
        if (cal.name != null && !cal.name.isEmpty()) {
            result.calibrations.add(cal);
        }
        for (ParsedCalibration child : cal.children) {
            collectCalibrations(child, result);
        }
    }

    static class ParseResult {
        ParsedCalibration calibration;
        int endIndex;
        ParseResult(ParsedCalibration cal, int end) {
            this.calibration = cal;
            this.endIndex = end;
        }
    }

    private static ParseResult parseSubtree(String newick, int start) {
        int i = start;
        while (i < newick.length() && Character.isWhitespace(newick.charAt(i))) i++;

        if (i >= newick.length()) {
            return new ParseResult(null, i);
        }

        if (newick.charAt(i) == '(') {
            return parseInternalNode(newick, i);
        } else {
            return parseLeaf(newick, i);
        }
    }

    private static ParseResult parseInternalNode(String newick, int start) {
        int i = start + 1;
        ParsedCalibration cal = new ParsedCalibration();

        while (i < newick.length()) {
            while (i < newick.length() && Character.isWhitespace(newick.charAt(i))) i++;
            if (i >= newick.length()) break;

            char c = newick.charAt(i);
            if (c == ')') {
                i++;
                break;
            } else if (c == ',') {
                i++;
                continue;
            } else {
                ParseResult child = parseSubtree(newick, i);
                if (child.calibration != null) {
                    cal.taxa.addAll(child.calibration.taxa);
                    if (child.calibration.name != null && !child.calibration.name.isEmpty()) {
                        cal.children.add(child.calibration);
                    }
                }
                i = child.endIndex;
            }
        }

        if (i < newick.length()) {
            Matcher m = ANNOTATION_PATTERN.matcher(newick.substring(i));
            if (m.lookingAt()) {
                applyAnnotation(cal, m.group(1));
                i += m.end();
            }
        }

        return new ParseResult(cal, i);
    }

    private static ParseResult parseLeaf(String newick, int start) {
        int i = start;
        StringBuilder name = new StringBuilder();

        while (i < newick.length()) {
            char c = newick.charAt(i);
            if (c == ',' || c == ')' || c == ':' || c == '[') break;
            name.append(c);
            i++;
        }

        String taxonName = name.toString().trim();
        ParsedCalibration cal = new ParsedCalibration();
        if (!taxonName.isEmpty()) {
            cal.taxa.add(taxonName);
        }

        if (i < newick.length() && newick.charAt(i) == '[') {
            Matcher m = ANNOTATION_PATTERN.matcher(newick.substring(i));
            if (m.lookingAt()) {
                applyAnnotation(cal, m.group(1));
                i += m.end();
            }
        }

        if (i < newick.length() && newick.charAt(i) == ':') {
            i++;
            while (i < newick.length()) {
                char c = newick.charAt(i);
                if (c == ',' || c == ')' || c == '[') break;
                i++;
            }
        }

        return new ParseResult(cal, i);
    }

    private static void applyAnnotation(ParsedCalibration cal, String annotation) {
        Matcher m = KEY_VALUE_PATTERN.matcher(annotation);
        while (m.find()) {
            String key = m.group(1).trim().toLowerCase();
            String value = m.group(2).trim();

            switch (key) {
                case "name":
                    cal.name = value;
                    break;
                case "dist":
                    cal.distType = value;
                    break;
                case "m":
                    try { cal.distM = Double.parseDouble(value); } catch (NumberFormatException ignored) {}
                    break;
                case "s":
                    try { cal.distS = Double.parseDouble(value); } catch (NumberFormatException ignored) {}
                    break;
                case "offset":
                    try { cal.distOffset = Double.parseDouble(value); } catch (NumberFormatException ignored) {}
                    break;
                case "monophyletic":
                    cal.monophyletic = Boolean.parseBoolean(value);
                    break;
            }
        }
    }
}
