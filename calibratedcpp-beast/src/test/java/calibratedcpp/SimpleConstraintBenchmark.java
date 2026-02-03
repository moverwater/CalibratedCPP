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
        System.out.printf("%-45s %6s %5s %10s %10s %10s%n",
                "File", "Taxa", "Cal", "p50(us)", "p95(us)", "p99(us)");
        System.out.println("-".repeat(95));

        for (Path newickFile : newickFiles) {
            try {
                // Skip files known to cause issues (overlapping calibrations, too large)
                String fileName = newickFile.getFileName().toString();
                if (fileName.contains("harrington") || fileName.contains("hara")) {
                    System.out.printf("%-45s %6s %5s %8s  SKIPPED (known issues)%n",
                            fileName.substring(0, Math.min(45, fileName.length())), "-", "-", "-");
                    continue;
                }

                BenchmarkResult result = runBenchmark(newickFile, iterations, warmup);
                System.out.printf("%-45s %6d %5d %10.1f %10.1f %10.1f%n",
                        fileName.substring(0, Math.min(45, fileName.length())),
                        result.numTaxa,
                        result.numCalibrations,
                        result.p50TimeUs,
                        result.p95TimeUs,
                        result.p99TimeUs);
            } catch (Exception e) {
                String msg = e.getMessage();
                if (msg != null && msg.length() > 60) msg = msg.substring(0, 60) + "...";
                System.out.printf("%-45s %6s %5s %8s  ERROR: %s%n",
                        newickFile.getFileName().toString().substring(0, Math.min(45, newickFile.getFileName().toString().length())),
                        "-", "-", "-", msg);
            }
        }
    }

    static class BenchmarkResult {
        int numTaxa;
        int numCalibrations;
        double meanTimeUs;
        double minTimeUs;
        double maxTimeUs;
        double p50TimeUs;
        double p95TimeUs;
        double p99TimeUs;
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

        // Generate random tree that respects calibration constraints (monophyletic)
        Tree tree = generateConstraintCompatibleTree(parsed.allTaxa, parsed.calibrations);

        // Debug: check monophyly
        // System.err.println("Generated tree for " + newickPath.getFileName());
        // System.err.println("  Num calibrations: " + calibrationClades.size());

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
        Arrays.sort(times);
        long sum = 0;
        for (long t : times) {
            sum += t;
        }

        result.meanTimeUs = (sum / (double) iterations) / 1000.0;
        result.minTimeUs = times[0] / 1000.0;
        result.maxTimeUs = times[iterations - 1] / 1000.0;
        result.p50TimeUs = times[iterations / 2] / 1000.0;
        result.p95TimeUs = times[(int)(iterations * 0.95)] / 1000.0;
        result.p99TimeUs = times[(int)(iterations * 0.99)] / 1000.0;

        return result;
    }

    /**
     * Generate a tree that respects calibration constraints:
     * 1. All calibration taxa are monophyletic
     * 2. Parent calibration MRCA is older than child calibration MRCA
     *
     * Strategy: Build tree directly from calibration hierarchy with proper heights.
     */
    private static Tree generateConstraintCompatibleTree(Set<String> taxaNames, List<ParsedCalibration> calibrations) throws Exception {
        // Build calibration hierarchy (find roots and children relationships)
        Map<ParsedCalibration, List<ParsedCalibration>> children = new HashMap<>();
        Map<ParsedCalibration, ParsedCalibration> parent = new HashMap<>();

        // Sort by size descending to process larger (parent) calibrations first
        List<ParsedCalibration> sorted = new ArrayList<>(calibrations);
        sorted.sort((a, b) -> b.taxa.size() - a.taxa.size());

        // Build parent-child relationships based on taxon set inclusion
        for (ParsedCalibration cal : sorted) {
            children.put(cal, new ArrayList<>());
        }

        for (int i = 0; i < sorted.size(); i++) {
            ParsedCalibration child = sorted.get(i);
            for (int j = 0; j < i; j++) {
                ParsedCalibration potentialParent = sorted.get(j);
                if (potentialParent.taxa.containsAll(child.taxa)) {
                    // Check if there's a closer parent
                    boolean hasCloserParent = false;
                    for (int k = j + 1; k < i; k++) {
                        ParsedCalibration middle = sorted.get(k);
                        if (middle.taxa.containsAll(child.taxa) && potentialParent.taxa.containsAll(middle.taxa)) {
                            hasCloserParent = true;
                            break;
                        }
                    }
                    if (!hasCloserParent) {
                        children.get(potentialParent).add(child);
                        parent.put(child, potentialParent);
                        break;
                    }
                }
            }
        }

        // Find root calibrations (no parent)
        List<ParsedCalibration> roots = new ArrayList<>();
        for (ParsedCalibration cal : calibrations) {
            if (!parent.containsKey(cal)) {
                roots.add(cal);
            }
        }

        // Assign depths (for height calculation)
        Map<ParsedCalibration, Integer> depth = new HashMap<>();
        assignDepths(roots, children, depth, 0);
        int maxDepth = depth.values().stream().mapToInt(Integer::intValue).max().orElse(0);

        // Build tree from calibration structure
        // Taxa not in any calibration
        Set<String> uncalibratedTaxa = new HashSet<>(taxaNames);
        for (ParsedCalibration cal : calibrations) {
            uncalibratedTaxa.removeAll(cal.taxa);
        }

        Random rand = new Random(42);
        double rootHeight = (maxDepth + 2) * 1.0;  // Root is oldest

        String newick = buildCalibrationTree(taxaNames, roots, children, depth, uncalibratedTaxa, rootHeight, rand);

        TreeParser treeParser = new TreeParser();
        treeParser.initByName("newick", newick + ";",
                "IsLabelledNewick", true,
                "adjustTipHeights", false);

        return treeParser;
    }

    private static void assignDepths(List<ParsedCalibration> nodes, Map<ParsedCalibration, List<ParsedCalibration>> children,
                                     Map<ParsedCalibration, Integer> depth, int currentDepth) {
        for (ParsedCalibration node : nodes) {
            depth.put(node, currentDepth);
            assignDepths(children.get(node), children, depth, currentDepth + 1);
        }
    }

    private static String buildCalibrationTree(Set<String> allTaxa, List<ParsedCalibration> roots,
                                                Map<ParsedCalibration, List<ParsedCalibration>> children,
                                                Map<ParsedCalibration, Integer> depth,
                                                Set<String> uncalibratedTaxa, double parentHeight, Random rand) {
        List<String> subtrees = new ArrayList<>();

        // Build subtrees for each root calibration
        for (ParsedCalibration root : roots) {
            String subtree = buildCalibrationSubtree(root, children, depth, parentHeight, rand);
            subtrees.add(subtree);
        }

        // Add uncalibrated taxa as individual leaves
        for (String taxon : uncalibratedTaxa) {
            subtrees.add(taxon + ":" + String.format("%.4f", parentHeight));
        }

        if (subtrees.isEmpty()) {
            return "";
        }
        if (subtrees.size() == 1) {
            return subtrees.get(0);
        }

        // Join all subtrees
        return "(" + String.join(",", subtrees) + "):" + String.format("%.4f", 0.1);
    }

    private static String buildCalibrationSubtree(ParsedCalibration cal, Map<ParsedCalibration, List<ParsedCalibration>> children,
                                                   Map<ParsedCalibration, Integer> depth, double parentHeight, Random rand) {
        int d = depth.get(cal);
        // Height decreases with depth (children are younger)
        double nodeHeight = parentHeight - 1.0 - rand.nextDouble() * 0.3;
        double branchLength = parentHeight - nodeHeight;

        List<ParsedCalibration> childCals = children.get(cal);

        // Find taxa directly in this calibration (not in any child calibration)
        Set<String> directTaxa = new HashSet<>(cal.taxa);
        for (ParsedCalibration child : childCals) {
            directTaxa.removeAll(child.taxa);
        }

        List<String> subtrees = new ArrayList<>();

        // Add child calibration subtrees
        for (ParsedCalibration child : childCals) {
            String childTree = buildCalibrationSubtree(child, children, depth, nodeHeight, rand);
            subtrees.add(childTree);
        }

        // Add direct taxa
        for (String taxon : directTaxa) {
            subtrees.add(taxon + ":" + String.format("%.4f", nodeHeight));
        }

        if (subtrees.size() == 1) {
            return subtrees.get(0);
        }

        // Build binary tree from subtrees
        return buildBinaryTree(subtrees, rand) + ":" + String.format("%.4f", branchLength);
    }

    private static String buildBinaryTree(List<String> subtrees, Random rand) {
        if (subtrees.size() == 1) {
            return subtrees.get(0);
        }
        if (subtrees.size() == 2) {
            return "(" + subtrees.get(0) + "," + subtrees.get(1) + ")";
        }

        // Split into two groups and recurse
        Collections.shuffle(subtrees, rand);
        int mid = subtrees.size() / 2;
        String left = buildBinaryTree(new ArrayList<>(subtrees.subList(0, mid)), rand);
        String right = buildBinaryTree(new ArrayList<>(subtrees.subList(mid, subtrees.size())), rand);

        return "(" + left + "," + right + ")";
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
