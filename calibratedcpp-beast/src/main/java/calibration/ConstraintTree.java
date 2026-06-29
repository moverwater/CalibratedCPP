package calibration;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibrationprior.CalibrationCladePrior;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parses clade definitions and optional upper/lower age bounds from a Newick string.
 * Distribution type is NOT stored here — it is a modelling choice on CalibrationPrior.
 *
 * @author Marcus Overwater
 */
@Description("Parses calibration clades (with optional upper/lower bounds) from a Newick string.")
public class ConstraintTree extends BEASTObject {

    public final Input<String> newickInput =
            new Input<>("newick",
                    "Constraint tree in Newick format. Internal nodes may carry "
                    + "[&lower=X,upper=Y] annotations. Empty string = no calibrations.",
                    "");

    private static final Pattern ANNOTATION_BLOCK = Pattern.compile("\\[&([^]]+)]");
    private static final Pattern KEY_VALUE         = Pattern.compile("([^=,]+)=([^,\\]]+)");

    private List<TaxonSet>              taxonSets;
    private List<CalibrationCladePrior> calibrationCladePriors;
    private Set<String>                 allTaxa;

    @Override
    public void initAndValidate() {
        String newick = newickInput.get().trim();
        if (newick.isEmpty() || newick.equals(";")) {
            taxonSets              = Collections.emptyList();
            calibrationCladePriors = Collections.emptyList();
            allTaxa                = Collections.emptySet();
            return;
        }

        ParsedNode root = parseNewick(newick);

        List<ParsedNode> calibrationNodes = new ArrayList<>();
        Set<String> taxa = new LinkedHashSet<>();
        collectNodes(root, calibrationNodes, taxa);

        this.allTaxa = Collections.unmodifiableSet(taxa);

        List<TaxonSet>              taxonSetList = new ArrayList<>();
        List<CalibrationCladePrior> priorList    = new ArrayList<>();

        for (ParsedNode node : calibrationNodes) {
            TaxonSet ts = buildTaxonSet(node);
            taxonSetList.add(ts);

            if (node.hasCalibrationBounds()) {
                CalibrationCladePrior prior = new CalibrationCladePrior();
                prior.initByName(
                        "taxa",     ts,
                        "lowerAge", new RealScalarParam<>(node.lower, NonNegativeReal.INSTANCE),
                        "upperAge", new RealScalarParam<>(node.upper, NonNegativeReal.INSTANCE)
                );
                if (node.name != null) prior.setID(node.name);
                priorList.add(prior);
            }
        }

        this.taxonSets              = Collections.unmodifiableList(taxonSetList);
        this.calibrationCladePriors = Collections.unmodifiableList(priorList);
    }

    // ── Public API ──────────────────────────────────────────────────────────────

    /** Taxon sets for each calibrated clade. */
    public List<TaxonSet> getTaxonSets() { return taxonSets; }

    /**
     * {@link CalibrationCladePrior} objects for clades that carry {@code lower} and {@code upper}
     * bounds. Pass to {@code CalibrationPrior} via its {@code calibration} input, or use
     * {@code CalibrationPrior}'s {@code constraintTree} input directly.
     */
    public List<CalibrationCladePrior> getCalibrationCladePriors() { return calibrationCladePriors; }

    /** All leaf taxon names in the tree. */
    public Set<String> getAllTaxa() { return allTaxa; }

    // ── Internal helpers ────────────────────────────────────────────────────────

    private static TaxonSet buildTaxonSet(ParsedNode node) {
        List<Taxon> taxonList = new ArrayList<>();
        for (String name : node.taxa) {
            Taxon t = new Taxon(name);
            t.setID(name);
            taxonList.add(t);
        }
        TaxonSet ts = new TaxonSet();
        ts.initByName("taxon", taxonList);
        if (node.name != null) ts.setID(node.name);
        return ts;
    }

    // ── Newick parsing ──────────────────────────────────────────────────────────

    private static class ParsedNode {
        String name = null;
        Set<String> taxa = new LinkedHashSet<>();
        double lower = Double.NaN, upper = Double.NaN;
        boolean virtualRoot = false, isLeaf = false;
        final List<ParsedNode> children = new ArrayList<>();

        boolean hasCalibrationBounds() {
            return !Double.isNaN(lower) && !Double.isNaN(upper);
        }
    }

    private static ParsedNode parseNewick(String input) {
        String s = stripComments(input).trim();
        if (s.endsWith(";")) s = s.substring(0, s.length() - 1).trim();
        if (s.isEmpty()) {
            ParsedNode empty = new ParsedNode();
            empty.virtualRoot = true;
            return empty;
        }
        int[] pos = {0};
        return parseSubtree(s, pos);
    }

    private static ParsedNode parseSubtree(String s, int[] pos) {
        skipWS(s, pos);
        return (pos[0] < s.length() && s.charAt(pos[0]) == '(')
                ? parseInternal(s, pos)
                : parseLeaf(s, pos);
    }

    private static ParsedNode parseInternal(String s, int[] pos) {
        pos[0]++;
        ParsedNode node = new ParsedNode();
        while (pos[0] < s.length()) {
            skipWS(s, pos);
            char c = s.charAt(pos[0]);
            if (c == ')') { pos[0]++; break; }
            if (c == ',') { pos[0]++; continue; }
            ParsedNode child = parseSubtree(s, pos);
            node.children.add(child);
            node.taxa.addAll(child.taxa);
        }
        skipLabel(s, pos);
        skipBranchLength(s, pos);
        applyAnnotation(node, readAnnotation(s, pos));
        return node;
    }

    private static ParsedNode parseLeaf(String s, int[] pos) {
        ParsedNode node = new ParsedNode();
        node.isLeaf = true;
        StringBuilder sb = new StringBuilder();
        while (pos[0] < s.length()) {
            char c = s.charAt(pos[0]);
            if (c == ',' || c == ')' || c == ':' || c == '[') break;
            sb.append(c); pos[0]++;
        }
        String name = sb.toString().trim();
        if (!name.isEmpty()) { node.taxa.add(name); node.name = name; }
        applyAnnotation(node, readAnnotation(s, pos));
        skipBranchLength(s, pos);
        return node;
    }

    private static void skipLabel(String s, int[] pos) {
        skipWS(s, pos);
        while (pos[0] < s.length()) {
            char c = s.charAt(pos[0]);
            if (c == ':' || c == '[' || c == ',' || c == ')' || c == ';') break;
            pos[0]++;
        }
    }

    private static void skipBranchLength(String s, int[] pos) {
        skipWS(s, pos);
        if (pos[0] < s.length() && s.charAt(pos[0]) == ':') {
            pos[0]++;
            while (pos[0] < s.length()) {
                char c = s.charAt(pos[0]);
                if (c == ',' || c == ')' || c == '[' || c == ';') break;
                pos[0]++;
            }
        }
    }

    private static void skipWS(String s, int[] pos) {
        while (pos[0] < s.length() && Character.isWhitespace(s.charAt(pos[0]))) pos[0]++;
    }

    private static String readAnnotation(String s, int[] pos) {
        skipWS(s, pos);
        if (pos[0] >= s.length() || s.charAt(pos[0]) != '[') return null;
        Matcher m = ANNOTATION_BLOCK.matcher(s).region(pos[0], s.length());
        if (m.lookingAt()) { pos[0] = m.end(); return m.group(1); }
        return null;
    }

    private static void applyAnnotation(ParsedNode node, String annotation) {
        if (annotation == null) return;
        Matcher m = KEY_VALUE.matcher(annotation);
        while (m.find()) {
            String key   = m.group(1).trim().toLowerCase(Locale.ROOT);
            String value = m.group(2).trim();
            switch (key) {
                case "name"               -> node.name = value;
                case "lower", "lowerage"  -> { try { node.lower = Double.parseDouble(value); } catch (NumberFormatException ignored) {} }
                case "upper", "upperage"  -> { try { node.upper = Double.parseDouble(value); } catch (NumberFormatException ignored) {} }
                case "virtualroot"        -> node.virtualRoot = Boolean.parseBoolean(value);
                // distribution-type keys are intentionally ignored — distributions are a
                // modelling choice on CalibrationPrior, not data stored in the constraint tree
            }
        }
    }

    private static void collectNodes(ParsedNode node, List<ParsedNode> out, Set<String> allTaxa) {
        allTaxa.addAll(node.taxa);
        for (ParsedNode child : node.children) collectNodes(child, out, allTaxa);
        if (!node.isLeaf && !node.virtualRoot && (node.name != null || node.hasCalibrationBounds())) {
            out.add(node);
        }
    }

    private static String stripComments(String input) {
        StringBuilder sb = new StringBuilder();
        for (String line : input.split("\n", -1)) {
            if (!line.trim().startsWith("#")) sb.append(line.trim());
        }
        return sb.toString();
    }
}
