package calibratedcpp.beauti;

import java.util.*;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.CompoundDistribution;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Beta;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.InverseGamma;
import beast.base.spec.inference.distribution.Laplace;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.Normal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedBirthDeathSkylineModel;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.text.FontWeight;

public class CalibrationPriorInputEditor extends InputEditor.Base {

    record DistDef(String name, List<ParamDef> params) {}
    record ParamDef(String key, String label, double defaultVal) {}

    static final List<DistDef> DISTS = List.of(
        new DistDef("Log-Normal", List.of(
            new ParamDef("M",      "M (log)", 1.0),
            new ParamDef("S",      "S (log)", 0.5),
            new ParamDef("offset", "Offset",  0.0))),
        new DistDef("Normal", List.of(
            new ParamDef("mean",   "Mean",   5.0),
            new ParamDef("sigma",  "Sigma",  1.0),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Exponential", List.of(
            new ParamDef("mean",   "Mean",   2.0),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Gamma", List.of(
            new ParamDef("alpha",  "Shape",  2.0),
            new ParamDef("beta",   "Rate",   0.5),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Inverse Gamma", List.of(
            new ParamDef("alpha",  "Alpha",  3.0),
            new ParamDef("beta",   "Beta",   2.0),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Beta", List.of(
            new ParamDef("alpha",  "Alpha",  1.0),
            new ParamDef("beta",   "Beta",   1.0))),
        new DistDef("Uniform", List.of(
            new ParamDef("lower",  "Lower",  0.0),
            new ParamDef("upper",  "Upper", 10.0))),
        new DistDef("Laplace", List.of(
            new ParamDef("mu",     "Mean",   5.0),
            new ParamDef("scale",  "Scale",  1.0),
            new ParamDef("offset", "Offset", 0.0)))
    );

    public CalibrationPriorInputEditor(BeautiDoc doc) { super(doc); }
    public CalibrationPriorInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibrationPrior.class; }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int listItemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, listItemNr, isExpandOption, addButtons);
        // Keep the label added by Base (matches other prior rows), remove only the editor widget
        pane.getChildren().removeIf(n -> !(n instanceof Label));
        VBox mainBox = new VBox(6);
        mainBox.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(mainBox, Priority.ALWAYS);
        pane.getChildren().add(mainBox);
        rebuild((CalibrationPrior) beastObject, mainBox);
    }

    private void rebuild(CalibrationPrior cp, VBox mainBox) {
        mainBox.getChildren().clear();
        String partition = partitionOf(cp);
        CalibratedBirthDeathSkylineModel cpp = getCPP(partition);

        if (cpp == null || cpp.calibrationsInput.get().isEmpty()) {
            mainBox.getChildren().add(new Label("No calibration clades defined. Use 'Manage calibrations' to add them, then click Refresh."));
            return;
        }

        List<TaxonSet> clades = new ArrayList<>(cpp.calibrationsInput.get());
        boolean mrcaMode = hasMRCAPriors(partition);

        ToggleGroup tg = new ToggleGroup();
        RadioButton calRb  = new RadioButton("CalibrationPrior");
        RadioButton mrcaRb = new RadioButton("MRCAPriors");
        calRb.setToggleGroup(tg);
        mrcaRb.setToggleGroup(tg);
        calRb.setSelected(!mrcaMode);
        mrcaRb.setSelected(mrcaMode);
        HBox modeRow = new HBox(12, new Label("Calibration Prior:"), calRb, mrcaRb);
        modeRow.setPadding(new Insets(4, 0, 8, 0));
        modeRow.setAlignment(Pos.CENTER_LEFT);
        mainBox.getChildren().add(modeRow);

        // Ensure cp.cladesInput and the prior CompoundDistribution stay consistent
        // on every render, not just when the user explicitly clicks a radio button.
        if (!mrcaMode) syncCalibrationPriorClades(cp, clades, partition);

        VBox content = FXUtils.newVBox();
        content.setSpacing(4);
        mainBox.getChildren().add(content);
        renderRows(content, cp, clades, partition, !mrcaMode);

        calRb.setOnAction(e -> {
            switchToCalibrationPrior(cp, clades, partition);
            renderRows(content, cp, clades, partition, true);
        });
        mrcaRb.setOnAction(e -> {
            switchToMRCAPriors(cp, clades, partition);
            renderRows(content, cp, clades, partition, false);
        });
    }

    // ── Row rendering ─────────────────────────────────────────────────────────────

    private void renderRows(VBox box, CalibrationPrior cp, List<TaxonSet> clades,
                             String partition, boolean calMode) {
        box.getChildren().clear();

        HBox hdr = new HBox(8);
        hdr.getChildren().add(fixedLabel("Clade", 160, true));
        if (calMode) {
            hdr.getChildren().addAll(fixedLabel("Lower bound", 100, true),
                                     fixedLabel("Upper bound", 100, true));
        } else {
            hdr.getChildren().addAll(fixedLabel("Distribution", 130, true),
                                     fixedLabel("Parameters", 300, true),
                                     fixedLabel("Monophyletic", 90, true));
        }
        box.getChildren().add(hdr);
        box.getChildren().add(new Separator());

        for (TaxonSet ts : clades) {
            HBox row = new HBox(8);
            row.setPadding(new Insets(3, 0, 3, 0));
            row.setAlignment(Pos.CENTER_LEFT);
            row.getChildren().add(fixedLabel(labelOf(ts), 160, false));
            if (calMode) buildCalRow(row, cp, ts, partition);
            else         buildMRCARow(row, cp, ts, partition);
            box.getChildren().add(row);
        }
    }

    private void buildCalRow(HBox row, CalibrationPrior cp, TaxonSet ts, String partition) {
        CalibrationCladePrior ccp = findCladePrior(ts, partition);
        TextField loTf = numField(ccp != null ? fmt(ccp.getLower()) : "");
        TextField hiTf = numField(ccp != null ? fmt(ccp.getUpper()) : "");
        loTf.setPrefWidth(100); loTf.setPromptText("optional");
        hiTf.setPrefWidth(100); hiTf.setPromptText("optional");
        Runnable apply = () -> applyCladePrior(cp, ts, partition, loTf.getText(), hiTf.getText());
        loTf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
        hiTf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
        row.getChildren().addAll(loTf, hiTf);
    }

    private void buildMRCARow(HBox row, CalibrationPrior cp, TaxonSet ts, String partition) {
        MRCAPrior mrca = findMRCAPrior(ts, partition);
        DistDef curDef = findDef(detectDistName(mrca));

        ChoiceBox<String> distCb = new ChoiceBox<>();
        DISTS.forEach(d -> distCb.getItems().add(d.name()));
        distCb.setValue(curDef.name());
        distCb.setPrefWidth(130);

        HBox paramsBox = new HBox(6);
        paramsBox.setAlignment(Pos.CENTER_LEFT);
        renderParamFields(paramsBox, curDef, mrca);

        CheckBox monoCb = new CheckBox();
        monoCb.setSelected(mrca == null || Boolean.TRUE.equals(mrca.getInputValue("monophyletic")));

        Runnable apply = () -> applyMRCAPrior(cp, ts, partition, findDef(distCb.getValue()),
                                               readParams(paramsBox), monoCb.isSelected());
        attachParamListeners(paramsBox, apply);
        monoCb.setOnAction(e -> apply.run());

        distCb.setOnAction(e -> {
            DistDef nd = findDef(distCb.getValue());
            renderParamFields(paramsBox, nd, null);
            attachParamListeners(paramsBox, apply);
            applyMRCAPrior(cp, ts, partition, nd, defaultParams(nd), monoCb.isSelected());
        });

        row.getChildren().addAll(distCb, paramsBox, monoCb);
    }

    private void renderParamFields(HBox paramsBox, DistDef def, MRCAPrior existingMrca) {
        paramsBox.getChildren().clear();
        for (ParamDef pd : def.params()) {
            double val = getDistParam(existingMrca, def.name(), pd.key(), pd.defaultVal());
            Label lbl = new Label(pd.label() + ":");
            lbl.setStyle("-fx-font-size:11px");
            TextField tf = numField(fmt(val));
            tf.setPrefWidth(72);
            tf.setUserData(pd.key());
            paramsBox.getChildren().addAll(lbl, tf);
        }
    }

    private void attachParamListeners(HBox paramsBox, Runnable apply) {
        for (var node : paramsBox.getChildren()) {
            if (node instanceof TextField tf) {
                tf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
            }
        }
    }

    private Map<String, Double> readParams(HBox paramsBox) {
        Map<String, Double> result = new LinkedHashMap<>();
        for (var node : paramsBox.getChildren()) {
            if (node instanceof TextField tf && tf.getUserData() instanceof String key) {
                try { result.put(key, Double.parseDouble(tf.getText().trim())); } catch (NumberFormatException ignored) {}
            }
        }
        return result;
    }

    private static Map<String, Double> defaultParams(DistDef def) {
        Map<String, Double> m = new LinkedHashMap<>();
        def.params().forEach(pd -> m.put(pd.key(), pd.defaultVal()));
        return m;
    }

    // ── Mode switching ────────────────────────────────────────────────────────────

    private void switchToCalibrationPrior(CalibrationPrior cp, List<TaxonSet> clades, String partition) {
        removeMRCAPriors(cp, partition);
        syncCalibrationPriorClades(cp, clades, partition);
    }

    private void syncCalibrationPriorClades(CalibrationPrior cp, List<TaxonSet> clades, String partition) {
        cp.cladesInput.get().clear();
        for (TaxonSet ts : clades) {
            CalibrationCladePrior ccp = findCladePrior(ts, partition);
            if (ccp != null) cp.cladesInput.get().add(ccp);
        }
        if (doc.pluginmap.get("prior") instanceof CompoundDistribution cd
                && !cd.pDistributions.get().contains(cp))
            cd.pDistributions.get().add(cp);
    }

    private void switchToMRCAPriors(CalibrationPrior cp, List<TaxonSet> clades, String partition) {
        cp.cladesInput.get().clear();
        cp.mrcaPriorsInput.get().clear();
        DistDef uniform = findDef("Uniform");
        for (TaxonSet ts : clades) {
            Map<String, Double> params = defaultParams(uniform);
            CalibrationCladePrior ccp = findCladePrior(ts, partition);
            if (ccp != null) {
                params.put("lower", ccp.getLower());
                params.put("upper", ccp.getUpper());
            }
            applyMRCAPrior(cp, ts, partition, uniform, params, true);
        }
    }

    private void removeMRCAPriors(CalibrationPrior cp, String partition) {
        cp.mrcaPriorsInput.get().clear();
        String pfx = "MRCAPrior.", sfx = "." + partition;
        List<String> ids = doc.pluginmap.keySet().stream()
            .filter(id -> id.startsWith(pfx) && id.endsWith(sfx))
            .toList();
        ids.forEach(doc.pluginmap::remove);
        List<String> distIds = doc.pluginmap.keySet().stream()
            .filter(id -> id.startsWith("MRCAPriorDist.") && id.endsWith(sfx))
            .toList();
        distIds.forEach(doc.pluginmap::remove);
    }

    // ── Apply actions ─────────────────────────────────────────────────────────────

    private void applyCladePrior(CalibrationPrior cp, TaxonSet ts, String partition,
                                  String loStr, String hiStr) {
        String ccpId = "CalibrationCladePrior." + labelOf(ts) + "." + partition;
        CalibrationCladePrior existing = (doc.pluginmap.get(ccpId) instanceof CalibrationCladePrior c) ? c : null;
        cp.cladesInput.get().remove(existing);

        loStr = loStr.trim(); hiStr = hiStr.trim();
        if (loStr.isEmpty() || hiStr.isEmpty()) {
            if (existing != null) doc.pluginmap.remove(ccpId);
            return;
        }
        try {
            double lo = Double.parseDouble(loStr);
            double hi = Double.parseDouble(hiStr);
            CalibrationCladePrior ccp;
            if (existing != null) {
                ccp = existing;
                ccp.lowerAgeInput.setValue(new RealScalarParam<>(lo, NonNegativeReal.INSTANCE), ccp);
                ccp.upperAgeInput.setValue(new RealScalarParam<>(hi, NonNegativeReal.INSTANCE), ccp);
                try { ccp.initAndValidate(); } catch (Exception ignored) {}
            } else {
                ccp = new CalibrationCladePrior();
                ccp.initByName("taxa", ts,
                    "lowerAge", new RealScalarParam<>(lo, NonNegativeReal.INSTANCE),
                    "upperAge", new RealScalarParam<>(hi, NonNegativeReal.INSTANCE));
                ccp.setID(ccpId);
                doc.addPlugin(ccp);
            }
            if (!cp.cladesInput.get().contains(ccp))
                cp.cladesInput.get().add(ccp);
        } catch (Exception ignored) {}
    }

    private void applyMRCAPrior(CalibrationPrior cp, TaxonSet ts, String partition, DistDef def,
                                  Map<String, Double> params, boolean monophyletic) {
        CalibratedBirthDeathSkylineModel cpp = getCPP(partition);
        if (cpp == null) return;

        String mrcaId = "MRCAPrior." + labelOf(ts) + "." + partition;
        String distId = "MRCAPriorDist." + labelOf(ts) + "." + partition;

        // Remove old objects from pluginmap first so addPlugin won't see a conflict
        if (doc.pluginmap.get(mrcaId) instanceof MRCAPrior old)
            cp.mrcaPriorsInput.get().remove(old);
        doc.pluginmap.remove(mrcaId);
        doc.pluginmap.remove(distId);

        ScalarDistribution dist = buildDistribution(def, params);
        if (dist != null) { dist.setID(distId); doc.addPlugin(dist); }

        MRCAPrior mrca = new MRCAPrior();
        mrca.setInputValue("tree",        cpp.treeInput.get());
        mrca.setInputValue("taxonset",    ts);
        mrca.setInputValue("monophyletic", monophyletic);
        if (dist != null) mrca.setInputValue("distr", dist);
        try { mrca.initAndValidate(); } catch (Exception e) { e.printStackTrace(); }
        mrca.setID(mrcaId);
        doc.addPlugin(mrca);
        // Add to CalibrationPrior's internal list, NOT to the top-level prior compound,
        // so MRCA priors don't appear as separate rows in BEAUti's Priors panel.
        if (!cp.mrcaPriorsInput.get().contains(mrca))
            cp.mrcaPriorsInput.get().add(mrca);
    }

    // ── Distribution builder ──────────────────────────────────────────────────────

    private ScalarDistribution buildDistribution(DistDef def, Map<String, Double> params) {
        try {
            return switch (def.name()) {
                case "Log-Normal" -> {
                    LogNormal d = new LogNormal();
                    d.setInputValue("M",      new RealScalarParam<>(params.getOrDefault("M",      1.0), Real.INSTANCE));
                    d.setInputValue("S",      new RealScalarParam<>(params.getOrDefault("S",      0.5), PositiveReal.INSTANCE));
                    d.setInputValue("offset", new RealScalarParam<>(params.getOrDefault("offset", 0.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Normal" -> {
                    Normal d = new Normal();
                    d.setInputValue("mean",   new RealScalarParam<>(params.getOrDefault("mean",   5.0), Real.INSTANCE));
                    d.setInputValue("sigma",  new RealScalarParam<>(params.getOrDefault("sigma",  1.0), PositiveReal.INSTANCE));
                    d.setInputValue("offset", new RealScalarParam<>(params.getOrDefault("offset", 0.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Exponential" -> {
                    Exponential d = new Exponential();
                    d.setInputValue("mean",   new RealScalarParam<>(params.getOrDefault("mean",   2.0), PositiveReal.INSTANCE));
                    d.setInputValue("offset", new RealScalarParam<>(params.getOrDefault("offset", 0.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Gamma" -> {
                    Gamma d = new Gamma();
                    d.setInputValue("alpha",  new RealScalarParam<>(params.getOrDefault("alpha",  2.0), PositiveReal.INSTANCE));
                    d.setInputValue("beta",   new RealScalarParam<>(params.getOrDefault("beta",   0.5), PositiveReal.INSTANCE));
                    d.setInputValue("offset", new RealScalarParam<>(params.getOrDefault("offset", 0.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Inverse Gamma" -> {
                    InverseGamma d = new InverseGamma();
                    d.setInputValue("alpha",  new RealScalarParam<>(params.getOrDefault("alpha",  3.0), PositiveReal.INSTANCE));
                    d.setInputValue("beta",   new RealScalarParam<>(params.getOrDefault("beta",   2.0), PositiveReal.INSTANCE));
                    d.setInputValue("offset", new RealScalarParam<>(params.getOrDefault("offset", 0.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Beta" -> {
                    Beta d = new Beta();
                    d.setInputValue("alpha", new RealScalarParam<>(params.getOrDefault("alpha", 1.0), PositiveReal.INSTANCE));
                    d.setInputValue("beta",  new RealScalarParam<>(params.getOrDefault("beta",  1.0), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Uniform" -> {
                    Uniform d = new Uniform();
                    d.setInputValue("lower", new RealScalarParam<>(params.getOrDefault("lower",  0.0), Real.INSTANCE));
                    d.setInputValue("upper", new RealScalarParam<>(params.getOrDefault("upper", 10.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Laplace" -> {
                    Laplace d = new Laplace();
                    d.setInputValue("mu",    new RealScalarParam<>(params.getOrDefault("mu",    5.0), Real.INSTANCE));
                    d.setInputValue("scale", new RealScalarParam<>(params.getOrDefault("scale", 1.0), PositiveReal.INSTANCE));
                    d.setInputValue("offset",new RealScalarParam<>(params.getOrDefault("offset",0.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                default -> null;
            };
        } catch (Exception e) { e.printStackTrace(); return null; }
    }

    // ── Lookups ───────────────────────────────────────────────────────────────────

    private CalibrationCladePrior findCladePrior(TaxonSet ts, String partition) {
        String id = "CalibrationCladePrior." + labelOf(ts) + "." + partition;
        return (doc.pluginmap.get(id) instanceof CalibrationCladePrior c) ? c : null;
    }

    private MRCAPrior findMRCAPrior(TaxonSet ts, String partition) {
        String id = "MRCAPrior." + labelOf(ts) + "." + partition;
        return (doc.pluginmap.get(id) instanceof MRCAPrior m) ? m : null;
    }

    private boolean hasMRCAPriors(String partition) {
        String sfx = "." + partition;
        return doc.pluginmap.keySet().stream().anyMatch(id -> id.startsWith("MRCAPrior.") && id.endsWith(sfx));
    }

    private String detectDistName(MRCAPrior mrca) {
        if (mrca == null) return "Uniform";
        Object d = mrca.getInputValue("distr");
        if (d instanceof LogNormal)    return "Log-Normal";
        if (d instanceof Normal)       return "Normal";
        if (d instanceof Exponential)  return "Exponential";
        if (d instanceof Gamma)        return "Gamma";
        if (d instanceof InverseGamma) return "Inverse Gamma";
        if (d instanceof Beta)         return "Beta";
        if (d instanceof Uniform)      return "Uniform";
        if (d instanceof Laplace)      return "Laplace";
        return "Uniform";
    }

    private double getDistParam(MRCAPrior mrca, String distName, String key, double defaultVal) {
        if (mrca == null) return defaultVal;
        Object distr = mrca.getInputValue("distr");
        if (!(distr instanceof BEASTInterface bi)) return defaultVal;
        try {
            Object v = bi.getInputValue(key);
            if (v instanceof Number n)            return n.doubleValue();
            if (v instanceof RealScalarParam<?> r) return (double) r.get(0);
        } catch (Exception ignored) {}
        return defaultVal;
    }

    private DistDef findDef(String name) {
        return DISTS.stream().filter(d -> d.name().equals(name)).findFirst()
               .orElse(DISTS.stream().filter(d -> d.name().equals("Uniform")).findFirst().orElse(DISTS.get(0)));
    }

    private static String labelOf(TaxonSet ts) {
        String id = ts.getID();
        if (id == null) return "unknown";
        // ID format: TaxonSet.{label}.{partition}
        int first = id.indexOf('.');
        int second = first >= 0 ? id.indexOf('.', first + 1) : -1;
        return second >= 0 ? id.substring(first + 1, second) : id;
    }

    private static String partitionOf(CalibrationPrior cp) {
        return cp.getID().replaceFirst("^CalibrationPrior\\.", "");
    }

    private CalibratedBirthDeathSkylineModel getCPP(String partition) {
        return (doc.pluginmap.get("CalibratedCPP." + partition) instanceof CalibratedBirthDeathSkylineModel m) ? m : null;
    }

    // ── UI helpers ────────────────────────────────────────────────────────────────

    private static Label fixedLabel(String text, double width, boolean bold) {
        Label lbl = new Label(text);
        lbl.setMinWidth(width);
        lbl.setPrefWidth(width);
        if (bold) lbl.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        else      lbl.setStyle("-fx-font-size: 12px;");
        return lbl;
    }

    private static TextField numField(String text) {
        TextField tf = new TextField(text);
        tf.setStyle("-fx-font-size: 12px;");
        return tf;
    }

    private static String fmt(double v) {
        if (Double.isNaN(v)) return "";
        return v == Math.floor(v) && !Double.isInfinite(v) ? String.valueOf((long) v) : String.valueOf(v);
    }
}
