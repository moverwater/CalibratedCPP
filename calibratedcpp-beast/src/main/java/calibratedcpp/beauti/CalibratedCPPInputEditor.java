package calibratedcpp.beauti;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.Tensor;
import beastfx.app.beauti.TreeDistributionInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedBirthDeathSkylineModel;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.SkylineParameter;
import calibration.ConstraintTree;
import calibrationprior.CalibrationCladePrior;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.stage.Modality;
import javafx.stage.Stage;

public class CalibratedCPPInputEditor extends TreeDistributionInputEditor {

    // ── Calibration state ─────────────────────────────────────────────────────────

    static class CalibrationEntry {
        String label;
        final String color;
        final List<String> directTaxa = new ArrayList<>();
        final List<CalibrationEntry> directChildCals = new ArrayList<>();
        Double lower, upper;
        // checkbox references — populated when the card is built, used for in-place updates
        final Map<String, CheckBox> taxaCbMap = new HashMap<>();
        final Map<CalibrationEntry, CheckBox> childCalCbMap = new HashMap<>();
        CalibrationEntry(String label, String color) { this.label = label; this.color = color; }
    }

    private static final String[] COLORS = {
        "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
        "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac"
    };

    private List<String> taxa = Arrays.asList("Apple", "Banana", "Cherry");
    private final List<CalibrationEntry> calibrationEntries = new ArrayList<>();
    private int colorIndex = 0;

    private VBox skylineEditorsBox;

    // ── Parameterization tables ───────────────────────────────────────────────────

    private static final String[][] PARAM_PAIRS = {
        {"Diversification rate + Turnover",            "diversificationRate", "turnover"},
        {"Birth rate + Death rate",                    "birthRate",           "deathRate"},
        {"Birth rate + Turnover",                      "birthRate",           "turnover"},
        {"Birth rate + Reproductive number",           "birthRate",           "reproductiveNumber"},
        {"Diversification rate + Death rate",          "diversificationRate", "deathRate"},
        {"Diversification rate + Reproductive number", "diversificationRate", "reproductiveNumber"},
        {"Death rate + Turnover",                      "deathRate",           "turnover"},
        {"Death rate + Reproductive number",           "deathRate",           "reproductiveNumber"},
        {"Birth rate + Diversification rate",          "birthRate",           "diversificationRate"},
    };

    private static final String[] ALL_RATE_NAMES =
        {"birthRate", "deathRate", "diversificationRate", "reproductiveNumber", "turnover"};

    public CalibratedCPPInputEditor(BeautiDoc doc) { super(doc); }
    public CalibratedCPPInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibratedCoalescentPointProcess.class; }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int listItemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, listItemNr, isExpandOption, addButtons);

        CalibratedBirthDeathSkylineModel cpp = (CalibratedBirthDeathSkylineModel) m_beastObject;
        try { taxa = cpp.treeInput.get().getTaxonset().asStringList(); } catch (Exception ignored) {}

        if (calibrationEntries.isEmpty()) initCalibrationsFromModel(cpp);

        Button calibrationButton = new Button("Manage calibrations");
        calibrationButton.setOnAction(e -> openPopup());
        pane.getChildren().add(calibrationButton);

        HBox paramRow = FXUtils.newHBox();
        paramRow.setSpacing(8);
        paramRow.setPadding(new Insets(8, 0, 4, 0));
        paramRow.getChildren().add(new Label("Parameterization:"));
        ChoiceBox<String> paramChoice = new ChoiceBox<>();
        for (String[] pair : PARAM_PAIRS) paramChoice.getItems().add(pair[0]);
        paramChoice.getSelectionModel().select(detectCurrentPairIndex(cpp));
        paramRow.getChildren().add(paramChoice);
        pane.getChildren().add(paramRow);

        skylineEditorsBox = FXUtils.newVBox();
        skylineEditorsBox.setSpacing(10);
        pane.getChildren().add(skylineEditorsBox);
        refreshSkylineEditors(cpp);

        paramChoice.getSelectionModel().selectedIndexProperty().addListener((obs, oldIdx, newIdx) -> {
            if (newIdx.intValue() < 0) return;
            String[] pair = PARAM_PAIRS[newIdx.intValue()];
            try { switchParameterization(cpp, pair[1], pair[2]); } catch (Exception ex) { ex.printStackTrace(); }
            refreshSkylineEditors(cpp);
        });
    }

    // ── Model initialization ──────────────────────────────────────────────────────

    private void initCalibrationsFromModel(CalibratedBirthDeathSkylineModel cpp) {
        String partition = partitionOf(cpp);
        List<CalibrationCladePrior> ccps = doc.pluginmap.entrySet().stream()
            .filter(e -> e.getKey().startsWith("CalibrationCladePrior.") && e.getKey().endsWith("." + partition))
            .filter(e -> e.getValue() instanceof CalibrationCladePrior)
            .map(e -> (CalibrationCladePrior) e.getValue())
            .toList();
        if (ccps.isEmpty()) return;

        Set<String> claimed = new LinkedHashSet<>();
        for (CalibrationCladePrior ccp : ccps) {
            String key = ccp.getID();
            String label = key.substring("CalibrationCladePrior.".length(),
                                         key.length() - ("." + partition).length());
            CalibrationEntry entry = new CalibrationEntry(label, COLORS[colorIndex % COLORS.length]);
            colorIndex++;
            entry.lower = ccp.getLower();
            entry.upper = ccp.getUpper();
            TaxonSet ts = ccp.getTaxa();
            if (ts != null) {
                for (Taxon t : ts.getTaxonSet()) {
                    String name = t.getID();
                    if (taxa.contains(name) && !claimed.contains(name)) {
                        entry.directTaxa.add(name);
                        claimed.add(name);
                    }
                }
            }
            calibrationEntries.add(entry);
        }
    }

    // ── Parameterization ─────────────────────────────────────────────────────────

    private int detectCurrentPairIndex(CalibratedBirthDeathSkylineModel cpp) {
        List<String> active = new ArrayList<>();
        for (String name : ALL_RATE_NAMES)
            if (getSkylineInput(cpp, name).get() != null) active.add(name);
        if (active.size() == 2) {
            for (int i = 0; i < PARAM_PAIRS.length; i++) {
                if (PARAM_PAIRS[i][1].equals(active.get(0)) && PARAM_PAIRS[i][2].equals(active.get(1))) return i;
                if (PARAM_PAIRS[i][1].equals(active.get(1)) && PARAM_PAIRS[i][2].equals(active.get(0))) return i;
            }
        }
        return 0;
    }

    private void switchParameterization(CalibratedBirthDeathSkylineModel cpp, String p1, String p2) {
        for (String name : ALL_RATE_NAMES) {
            Input<SkylineParameter> in = getSkylineInput(cpp, name);
            SkylineParameter old = in.get();
            if (old != null) setEstimated(old.valuesInput.get(), false);
            in.setValue(null, cpp);
        }
        String partition = partitionOf(cpp);
        activateSkylineParam(cpp, p1, partition);
        activateSkylineParam(cpp, p2, partition);
    }

    private void activateSkylineParam(CalibratedBirthDeathSkylineModel cpp, String paramName, String partition) {
        String paramId = paramName + "." + partition;
        String spId    = paramName + "SP." + partition;
        RealScalarParam<?> scalar;
        BEASTInterface existing = doc.pluginmap.get(paramId);
        if (existing instanceof RealScalarParam<?> rsp) {
            scalar = rsp;
        } else {
            scalar = new RealScalarParam<>(defaultValueFor(paramName), domainFor(paramName));
            scalar.setID(paramId);
            pluginPut(paramId, scalar);
        }
        setEstimated(scalar, true);
        SkylineParameter sp;
        BEASTInterface existingSP = doc.pluginmap.get(spId);
        if (existingSP instanceof SkylineParameter esp) {
            sp = esp;
        } else {
            sp = new SkylineParameter();
            sp.valuesInput.setValue(scalar, sp);
            sp.initAndValidate();
            sp.setID(spId);
            pluginPut(spId, sp);
        }
        getSkylineInput(cpp, paramName).setValue(sp, cpp);
    }

    /** Registers a plugin, using doc.addPlugin when the ID is new, or directly replacing a stale
     *  placeholder when the ID already exists (avoiding BEAUti's "could not add entry" dialog). */
    private void pluginPut(String id, BEASTInterface plugin) {
        if (doc.pluginmap.containsKey(id))
            doc.pluginmap.put(id, plugin);
        else
            doc.addPlugin(plugin);
    }

    @SuppressWarnings("unchecked")
    private Input<SkylineParameter> getSkylineInput(CalibratedBirthDeathSkylineModel cpp, String name) {
        return (Input<SkylineParameter>) switch (name) {
            case "birthRate"           -> cpp.birthRateInput;
            case "deathRate"           -> cpp.deathRateInput;
            case "diversificationRate" -> cpp.diversificationRateInput;
            case "reproductiveNumber"  -> cpp.reproductiveNumberInput;
            case "turnover"            -> cpp.turnoverInput;
            default -> throw new IllegalArgumentException("Unknown rate: " + name);
        };
    }

    private static String partitionOf(CalibratedBirthDeathSkylineModel cpp) {
        return cpp.getID().replaceFirst("^CalibratedCPP\\.", "");
    }

    private static void setEstimated(Tensor<?, ?> tensor, boolean estimated) {
        if (tensor instanceof StateNode sn) sn.isEstimatedInput.setValue(estimated, sn);
    }

    // ── Skyline editors ───────────────────────────────────────────────────────────

    private void refreshSkylineEditors(CalibratedBirthDeathSkylineModel cpp) {
        skylineEditorsBox.getChildren().clear();
        for (String name : ALL_RATE_NAMES) {
            SkylineParameter sp = getSkylineInput(cpp, name).get();
            if (sp != null) skylineEditorsBox.getChildren().add(buildSkylineEditor(sp, name));
        }
    }

    private VBox buildSkylineEditor(SkylineParameter sp, String paramName) {
        VBox box = FXUtils.newVBox();
        box.setSpacing(4);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");
        box.getChildren().add(new Label(displayNameOf(paramName)));

        int nChanges = sp.changeTimesInput.get() == null ? 0 : sp.changeTimesInput.get().size();

        HBox epochRow = FXUtils.newHBox();
        epochRow.setSpacing(6);
        epochRow.getChildren().add(new Label("Epochs:"));
        Spinner<Integer> epochSpinner = new Spinner<>(1, 100, nChanges + 1);
        epochSpinner.setEditable(true);
        epochSpinner.setPrefWidth(70);
        epochRow.getChildren().add(epochSpinner);
        box.getChildren().add(epochRow);

        HBox flagsRow = FXUtils.newHBox();
        flagsRow.setSpacing(12);
        CheckBox agesCb = new CheckBox("Times are ages");
        agesCb.setSelected(sp.timesAreAgesInput.get());
        CheckBox relCb = new CheckBox("Relative to process length (0–1)");
        relCb.setSelected(sp.timesAreRelativeInput.get());
        relCb.setTooltip(new Tooltip("Change times as fractions of process length — must be between 0 and 1."));
        Button distBtn = new Button("Distribute evenly");
        distBtn.setTooltip(new Tooltip("Space change times evenly (fractions if relative, else 1, 2, ...)."));
        flagsRow.getChildren().addAll(agesCb, relCb, distBtn);
        box.getChildren().add(flagsRow);

        HBox ctRow = FXUtils.newHBox();
        ctRow.setSpacing(4);
        VBox ctBox = FXUtils.newVBox();
        ctBox.setSpacing(2);
        ctBox.getChildren().addAll(new Label("Change times:"), ctRow);
        ctBox.setVisible(nChanges > 0);
        ctBox.setManaged(nChanges > 0);
        box.getChildren().add(ctBox);

        HBox valRow = FXUtils.newHBox();
        valRow.setSpacing(4);
        VBox valBox = FXUtils.newVBox();
        valBox.setSpacing(2);
        valBox.getChildren().addAll(new Label("Epoch values (root → present):"), valRow);
        box.getChildren().add(valBox);

        rebuildSkylineRows(sp, nChanges, ctRow, valRow, agesCb, relCb);

        epochSpinner.valueProperty().addListener((obs, o, nv) -> {
            int nc = nv - 1;
            ensureChangeTimes(sp, nc);
            ensureValues(sp, nv);
            ctBox.setVisible(nc > 0);
            ctBox.setManaged(nc > 0);
            rebuildSkylineRows(sp, nc, ctRow, valRow, agesCb, relCb);
        });
        agesCb.selectedProperty().addListener((obs, o, n) -> {
            sp.timesAreAgesInput.setValue(n, sp);
            rebuildSkylineRows(sp, epochSpinner.getValue() - 1, ctRow, valRow, agesCb, relCb);
        });
        relCb.selectedProperty().addListener((obs, o, n) -> {
            sp.timesAreRelativeInput.setValue(n, sp);
            int nc = epochSpinner.getValue() - 1;
            ctBox.setVisible(nc > 0);
            ctBox.setManaged(nc > 0);
            rebuildSkylineRows(sp, nc, ctRow, valRow, agesCb, relCb);
            setTimeFieldValues(ctRow, evenlySpaced(nc, n));
            sp.initAndValidate();
            sync();
        });
        distBtn.setOnAction(e -> {
            int nc = epochSpinner.getValue() - 1;
            rebuildSkylineRows(sp, nc, ctRow, valRow, agesCb, relCb);
            setTimeFieldValues(ctRow, evenlySpaced(nc, relCb.isSelected()));
            sp.initAndValidate();
            sync();
        });
        return box;
    }

    private void rebuildSkylineRows(SkylineParameter sp, int nChanges,
                                    HBox ctRow, HBox valRow, CheckBox agesCb, CheckBox relCb) {
        ctRow.getChildren().clear();
        valRow.getChildren().clear();
        boolean relative = relCb.isSelected();
        if (nChanges > 0) {
            RealVectorParam<?> ctParam = changeTimesParam(sp);
            double[] defaults = evenlySpaced(nChanges, relative);
            for (int i = 0; i < nChanges; i++) {
                double v = (ctParam != null && i < ctParam.size()) ? (double) ctParam.get(i) : defaults[i];
                ctRow.getChildren().add(new Label(agesCb.isSelected() ? ("Age " + (i + 1) + ":") : ("t" + (i + 1) + ":")));
                final int fi = i;
                TextField tf = compactField(v);
                tf.textProperty().addListener((obs, old, nw) -> {
                    try { setChangeTime(sp, fi, Double.parseDouble(nw)); } catch (NumberFormatException ignored) {}
                });
                ctRow.getChildren().add(tf);
            }
        }
        int nEpochs = nChanges + 1;
        Tensor<?, Double> vParam = sp.valuesInput.get();
        for (int i = 0; i < nEpochs; i++) {
            double v = (vParam != null && i < vParam.size()) ? (double) vParam.get(i) : 1.0;
            valRow.getChildren().add(new Label("E" + (i + 1) + ":"));
            final int fi = i;
            TextField tf = compactField(v);
            tf.textProperty().addListener((obs, old, nw) -> {
                try { setEpochValue(sp, fi, Double.parseDouble(nw)); } catch (NumberFormatException ignored) {}
            });
            valRow.getChildren().add(tf);
        }
    }

    private static TextField compactField(double value) {
        TextField tf = new TextField(String.valueOf(value));
        tf.setPrefWidth(70);
        tf.setPadding(new Insets(2));
        return tf;
    }

    @SuppressWarnings("unchecked")
    private void ensureChangeTimes(SkylineParameter sp, int nChanges) {
        if (nChanges == 0) { sp.changeTimesInput.setValue(null, sp); return; }
        boolean relative = sp.timesAreRelativeInput.get();
        RealVectorParam<NonNegativeReal> ct = changeTimesParam(sp);
        if (ct == null) {
            // The vector may have been disconnected (setValue null) but still live in pluginmap —
            // reuse it to avoid a duplicate-ID dialog on the next addPlugin call.
            String spId = sp.getID();
            String ctId = spId != null ? spId + ".changeTimes" : null;
            if (ctId != null && doc.pluginmap.get(ctId) instanceof RealVectorParam<?> v)
                ct = (RealVectorParam<NonNegativeReal>) v;
            else {
                ct = new RealVectorParam<>(evenlySpaced(nChanges, relative), NonNegativeReal.INSTANCE);
                if (ctId != null) ct.setID(ctId);
                doc.addPlugin(ct);
            }
            sp.changeTimesInput.setValue(ct, sp);
        }
        double[] old = toDoubles(ct);
        double[] defaults = evenlySpaced(nChanges, relative);
        double[] upd = new double[nChanges];
        for (int i = 0; i < nChanges; i++) upd[i] = i < old.length ? old[i] : defaults[i];
        resetVector(ct, upd);
    }

    private static void setTimeFieldValues(HBox ctRow, double[] values) {
        int vi = 0;
        for (javafx.scene.Node node : ctRow.getChildren()) {
            if (node instanceof TextField tf && vi < values.length)
                tf.setText(String.valueOf(values[vi++]));
        }
    }

    private static double[] evenlySpaced(int n, boolean relative) {
        double[] v = new double[n];
        for (int i = 0; i < n; i++) v[i] = relative ? (double)(i + 1) / (n + 1) : (double)(i + 1);
        return v;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private void ensureValues(SkylineParameter sp, int nEpochs) {
        Tensor<?, Double> existing = sp.valuesInput.get();
        Real dom = existing instanceof RealScalarParam<?> rsp ? rsp.domainTypeInput.get()
                 : existing instanceof RealVectorParam<?> vp  ? vp.domainTypeInput.get()
                 : Real.INSTANCE;
        if (existing instanceof RealVectorParam<?> vp) {
            double[] old = toDoubles(vp);
            double[] upd = new double[nEpochs];
            for (int i = 0; i < nEpochs; i++) upd[i] = i < old.length ? old[i] : (old.length > 0 ? old[old.length - 1] : 1.0);
            resetVector(vp, upd);
        } else {
            double scalar = (existing != null) ? (double) existing.get(0) : 1.0;
            double[] init = new double[nEpochs];
            for (int i = 0; i < nEpochs; i++) init[i] = scalar;
            String spId = sp.getID();
            String vpId = spId != null ? spId + ".values" : null;
            // Reuse an orphaned vector from pluginmap rather than adding a duplicate.
            @SuppressWarnings({"unchecked","rawtypes"})
            RealVectorParam<?> vp = (vpId != null && doc.pluginmap.get(vpId) instanceof RealVectorParam<?> v)
                    ? v : new RealVectorParam(init, dom);
            if (vp.getID() == null && vpId != null) vp.setID(vpId);
            if (!doc.pluginmap.containsKey(vp.getID())) doc.addPlugin(vp);
            sp.valuesInput.setValue(vp, sp);
        }
    }

    private void setChangeTime(SkylineParameter sp, int idx, double val) {
        RealVectorParam<?> ct = changeTimesParam(sp);
        if (ct == null) return;
        double[] vals = toDoubles(ct);
        if (idx < vals.length) { vals[idx] = val; resetVector(ct, vals); }
    }

    private void setEpochValue(SkylineParameter sp, int idx, double val) {
        if (!(sp.valuesInput.get() instanceof RealVectorParam<?> vp)) return;
        double[] vals = toDoubles(vp);
        if (idx < vals.length) { vals[idx] = val; resetVector(vp, vals); }
    }

    @SuppressWarnings("unchecked")
    private static RealVectorParam<NonNegativeReal> changeTimesParam(SkylineParameter sp) {
        Tensor<?, ?> t = sp.changeTimesInput.get();
        return (t instanceof RealVectorParam<?>) ? (RealVectorParam<NonNegativeReal>) t : null;
    }

    private static double[] toDoubles(RealVectorParam<?> vp) {
        double[] out = new double[vp.size()];
        for (int i = 0; i < out.length; i++) out[i] = (double) vp.get(i);
        return out;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static void resetVector(RealVectorParam<?> vp, double[] vals) {
        List<Double> list = new ArrayList<>(vals.length);
        for (double v : vals) list.add(v);
        ((RealVectorParam<Real>) vp).valuesInput.setValue(list, vp);
        vp.initAndValidate();
    }

    private static String displayNameOf(String paramName) {
        return switch (paramName) {
            case "birthRate"           -> "Birth rate (λ)";
            case "deathRate"           -> "Death rate (μ)";
            case "diversificationRate" -> "Diversification rate (λ − μ)";
            case "turnover"            -> "Turnover (μ/λ)";
            case "reproductiveNumber"  -> "Reproductive number (λ/μ)";
            default -> paramName;
        };
    }

    private static double defaultValueFor(String paramName) {
        return switch (paramName) {
            case "birthRate"           -> 0.5;
            case "deathRate"           -> 0.25;
            case "diversificationRate" -> 0.1;
            case "turnover"            -> 0.5;
            case "reproductiveNumber"  -> 2.0;
            default -> 1.0;
        };
    }

    private static Real domainFor(String paramName) {
        return switch (paramName) {
            case "birthRate", "deathRate", "turnover", "reproductiveNumber" -> PositiveReal.INSTANCE;
            default -> Real.INSTANCE;
        };
    }

    // ── Calibration popup ─────────────────────────────────────────────────────────

    private void openPopup() {
        Stage stage = new Stage();
        stage.initModality(Modality.APPLICATION_MODAL);
        stage.setTitle("Manage Calibrations");

        VBox cardBox = FXUtils.newVBox();
        cardBox.setSpacing(8);
        cardBox.setPadding(new Insets(10));
        ScrollPane scrollPane = new ScrollPane(cardBox);
        scrollPane.setFitToWidth(true);
        scrollPane.setStyle("-fx-background-color: #f5f5f5;");

        HBox topBar = new HBox(8);
        topBar.setPadding(new Insets(8, 10, 8, 10));
        topBar.setAlignment(Pos.CENTER_LEFT);
        topBar.setStyle("-fx-background-color: #ececec; -fx-border-color: #d0d0d0; -fx-border-width: 0 0 1 0;");
        Label title = new Label("Calibrations");
        title.setStyle("-fx-font-weight: bold; -fx-font-size: 13px;");
        HBox.setHgrow(title, Priority.ALWAYS);
        Button addBtn = new Button("+ Add");
        addBtn.setStyle("-fx-background-color: #27ae60; -fx-text-fill: white; -fx-font-weight: bold;");
        topBar.getChildren().addAll(title, addBtn);

        HBox bottomBar = new HBox(8);
        bottomBar.setPadding(new Insets(8, 10, 8, 10));
        bottomBar.setAlignment(Pos.CENTER_LEFT);
        bottomBar.setStyle("-fx-background-color: #ececec; -fx-border-color: #d0d0d0; -fx-border-width: 1 0 0 0;");
        Button importBtn = new Button("Import Newick…");
        Button exportBtn = new Button("Export Newick…");
        Region spacer = new Region();
        HBox.setHgrow(spacer, Priority.ALWAYS);
        Button saveBtn = new Button("Save & Close");
        saveBtn.setStyle("-fx-background-color: #2980b9; -fx-text-fill: white; -fx-font-weight: bold;");
        bottomBar.getChildren().addAll(importBtn, exportBtn, spacer, saveBtn);

        Runnable rebuild = () -> rebuildCalibrationList(cardBox, scrollPane);

        addBtn.setOnAction(e -> {
            String color = COLORS[colorIndex % COLORS.length];
            colorIndex++;
            calibrationEntries.add(new CalibrationEntry("Clade_" + (calibrationEntries.size() + 1), color));
            rebuild.run();
        });
        importBtn.setOnAction(e -> importNewick(rebuild));
        exportBtn.setOnAction(e -> exportNewick());
        saveBtn.setOnAction(e -> { saveCalibrations(); stage.close(); });

        rebuild.run();

        VBox root = new VBox(topBar, scrollPane, bottomBar);
        VBox.setVgrow(scrollPane, Priority.ALWAYS);
        stage.setScene(new Scene(root, 680, 600));
        stage.showAndWait();
    }

    private void rebuildCalibrationList(VBox cardBox, ScrollPane scrollPane) {
        double savedV = scrollPane.getVvalue();
        cardBox.getChildren().clear();
        for (CalibrationEntry entry : calibrationEntries)
            cardBox.getChildren().add(buildCalibrationCard(entry, cardBox, scrollPane));
        Platform.runLater(() -> scrollPane.setVvalue(savedV));
    }

    private VBox buildCalibrationCard(CalibrationEntry entry, VBox cardBox, ScrollPane scrollPane) {
        VBox card = FXUtils.newVBox();
        card.setSpacing(6);
        card.setPadding(new Insets(10));
        card.setStyle(
            "-fx-background-color: white;" +
            "-fx-border-color: " + entry.color + " #ddd #ddd #ddd;" +
            "-fx-border-width: 3 1 1 1;" +
            "-fx-border-radius: 4;" +
            "-fx-background-radius: 4;"
        );

        // Header: [dot] [name field] [× delete]
        HBox headerRow = FXUtils.newHBox();
        headerRow.setSpacing(8);
        headerRow.setAlignment(Pos.CENTER_LEFT);
        Region dot = new Region();
        dot.setMinSize(12, 12); dot.setMaxSize(12, 12);
        dot.setStyle("-fx-background-color: " + entry.color + "; -fx-background-radius: 6;");
        TextField nameTf = new TextField(entry.label);
        HBox.setHgrow(nameTf, Priority.ALWAYS);
        nameTf.textProperty().addListener((obs, o, n) ->
            entry.label = n.replaceAll("[^\\w]", "_").isEmpty() ? entry.label : n.replaceAll("[^\\w]", "_"));
        Button delBtn = new Button("×");
        delBtn.setStyle("-fx-background-color: #e74c3c; -fx-text-fill: white; -fx-padding: 2 7;");
        delBtn.setOnAction(e -> {
            for (CalibrationEntry other : calibrationEntries) other.directChildCals.remove(entry);
            calibrationEntries.remove(entry);
            rebuildCalibrationList(cardBox, scrollPane);
        });
        headerRow.getChildren().addAll(dot, nameTf, delBtn);

        // Bounds row
        HBox boundsRow = FXUtils.newHBox();
        boundsRow.setSpacing(8);
        boundsRow.setAlignment(Pos.CENTER_LEFT);
        boundsRow.getChildren().add(new Label("Lower:"));
        TextField lowerTf = boundsField(entry.lower);
        lowerTf.textProperty().addListener((obs, o, n) -> {
            try { entry.lower = n.trim().isEmpty() ? null : Double.parseDouble(n); }
            catch (NumberFormatException ignored) {}
        });
        boundsRow.getChildren().add(lowerTf);
        boundsRow.getChildren().add(new Label("Upper:"));
        TextField upperTf = boundsField(entry.upper);
        upperTf.textProperty().addListener((obs, o, n) -> {
            try { entry.upper = n.trim().isEmpty() ? null : Double.parseDouble(n); }
            catch (NumberFormatException ignored) {}
        });
        boundsRow.getChildren().add(upperTf);

        card.getChildren().addAll(headerRow, boundsRow, new Separator());

        // Taxa section
        entry.taxaCbMap.clear();
        Label taxaHeading = new Label("Direct taxa:");
        taxaHeading.setStyle("-fx-font-size: 11px; -fx-font-weight: bold; -fx-text-fill: #555;");
        VBox taxaBox = new VBox(2);
        taxaBox.setPadding(new Insets(4, 4, 4, 8));
        for (String taxon : taxa) {
            boolean ownedByOther = calibrationEntries.stream()
                .filter(e -> e != entry).anyMatch(e -> e.directTaxa.contains(taxon));
            CheckBox cb = new CheckBox(taxon);
            cb.setSelected(entry.directTaxa.contains(taxon));
            cb.setDisable(ownedByOther && !entry.directTaxa.contains(taxon));
            if (ownedByOther && !entry.directTaxa.contains(taxon))
                cb.setTooltip(new Tooltip("Assigned to " + calibrationEntries.stream()
                    .filter(e -> e != entry && e.directTaxa.contains(taxon))
                    .map(e -> e.label).findFirst().orElse("?")));
            cb.selectedProperty().addListener((obs, o, n) -> {
                if (n) {
                    for (CalibrationEntry other : calibrationEntries)
                        if (other != entry) other.directTaxa.remove(taxon);
                    entry.directTaxa.add(taxon);
                } else {
                    entry.directTaxa.remove(taxon);
                }
                updateCheckBoxStates();
            });
            entry.taxaCbMap.put(taxon, cb);
            taxaBox.getChildren().add(cb);
        }
        ScrollPane taxaSp = new ScrollPane(taxaBox);
        taxaSp.setFitToWidth(true);
        taxaSp.setPrefViewportHeight(Math.min(taxa.size() * 24 + 8, 130));
        taxaSp.setStyle("-fx-background-color: #f9f9f9; -fx-border-color: #eee;");
        card.getChildren().addAll(taxaHeading, taxaSp);

        // Child calibrations section (only when there are others)
        entry.childCalCbMap.clear();
        List<CalibrationEntry> others = calibrationEntries.stream().filter(e -> e != entry).toList();
        if (!others.isEmpty()) {
            Label calHeading = new Label("Child calibrations:");
            calHeading.setStyle("-fx-font-size: 11px; -fx-font-weight: bold; -fx-text-fill: #555;");
            VBox calBox = new VBox(2);
            calBox.setPadding(new Insets(4, 4, 4, 8));
            for (CalibrationEntry other : others) {
                boolean ownedByOther = calibrationEntries.stream()
                    .filter(e -> e != entry).anyMatch(e -> e.directChildCals.contains(other));
                boolean wouldCycle = isAncestor(other, entry);
                CheckBox cb = new CheckBox(other.label);
                cb.setSelected(entry.directChildCals.contains(other));
                cb.setDisable((ownedByOther || wouldCycle) && !entry.directChildCals.contains(other));
                if (wouldCycle)
                    cb.setTooltip(new Tooltip("Would create a cycle"));
                else if (ownedByOther)
                    cb.setTooltip(new Tooltip("Child of " + calibrationEntries.stream()
                        .filter(e -> e != entry && e.directChildCals.contains(other))
                        .map(e -> e.label).findFirst().orElse("?")));
                cb.selectedProperty().addListener((obs, o, n) -> {
                    if (n) {
                        for (CalibrationEntry e : calibrationEntries)
                            if (e != entry) e.directChildCals.remove(other);
                        entry.directChildCals.add(other);
                    } else {
                        entry.directChildCals.remove(other);
                    }
                    updateCheckBoxStates();
                });
                entry.childCalCbMap.put(other, cb);
                calBox.getChildren().add(cb);
            }
            card.getChildren().addAll(calHeading, calBox);
        }

        return card;
    }

    private static boolean isAncestor(CalibrationEntry potentialAncestor, CalibrationEntry ofEntry) {
        for (CalibrationEntry child : ofEntry.directChildCals) {
            if (child == potentialAncestor) return true;
            if (isAncestor(potentialAncestor, child)) return true;
        }
        return false;
    }

    /** Updates disabled/tooltip state on all existing checkboxes without rebuilding any cards. */
    private void updateCheckBoxStates() {
        for (CalibrationEntry entry : calibrationEntries) {
            for (String taxon : taxa) {
                CheckBox cb = entry.taxaCbMap.get(taxon);
                if (cb == null) continue;
                boolean owned = entry.directTaxa.contains(taxon);
                boolean ownedByOther = calibrationEntries.stream()
                    .filter(e -> e != entry).anyMatch(e -> e.directTaxa.contains(taxon));
                cb.setDisable(ownedByOther && !owned);
                if (ownedByOther && !owned)
                    cb.setTooltip(new Tooltip("Assigned to " + calibrationEntries.stream()
                        .filter(e -> e != entry && e.directTaxa.contains(taxon))
                        .map(e -> e.label).findFirst().orElse("?")));
                else
                    cb.setTooltip(null);
            }
            for (CalibrationEntry other : calibrationEntries) {
                if (other == entry) continue;
                CheckBox cb = entry.childCalCbMap.get(other);
                if (cb == null) continue;
                boolean owned = entry.directChildCals.contains(other);
                boolean ownedByOther = calibrationEntries.stream()
                    .filter(e -> e != entry).anyMatch(e -> e.directChildCals.contains(other));
                boolean wouldCycle = isAncestor(other, entry);
                cb.setDisable((ownedByOther || wouldCycle) && !owned);
                if (wouldCycle) cb.setTooltip(new Tooltip("Would create a cycle"));
                else if (ownedByOther)
                    cb.setTooltip(new Tooltip("Child of " + calibrationEntries.stream()
                        .filter(e -> e != entry && e.directChildCals.contains(other))
                        .map(e -> e.label).findFirst().orElse("?")));
                else cb.setTooltip(null);
            }
        }
    }

    private static TextField boundsField(Double value) {
        TextField tf = new TextField(value == null ? "" : String.valueOf(value));
        tf.setPrefWidth(90);
        tf.setPromptText("optional");
        return tf;
    }

    // ── Import / Export Newick ────────────────────────────────────────────────────

    private void importNewick(Runnable rebuildFn) {
        Stage dialog = new Stage();
        dialog.initModality(Modality.APPLICATION_MODAL);
        dialog.setTitle("Import Constraint Tree (Newick)");

        TextArea ta = new TextArea();
        ta.setPromptText("Paste Newick constraint tree here, e.g.:\n((A,B)[&name=Clade1,lower=5,upper=10],C)[&name=Root,lower=15,upper=30];");
        ta.setWrapText(true);
        ta.setPrefRowCount(6);

        Label hint = new Label("Annotations: [&name=Label,lower=X,upper=Y]  — name and bounds are optional.");
        hint.setStyle("-fx-font-size:11px; -fx-text-fill:#555;");

        Button okBtn = new Button("Import");
        Button cancelBtn = new Button("Cancel");
        cancelBtn.setOnAction(e -> dialog.close());
        okBtn.setOnAction(e -> {
            String newick = ta.getText().trim();
            if (newick.isEmpty()) { dialog.close(); return; }
            try {
                ConstraintTree ct = new ConstraintTree();
                ct.initByName("newick", newick);
                List<TaxonSet> allTs = ct.getTaxonSets();

                Map<TaxonSet, Set<String>> allTaxaOf = new HashMap<>();
                for (TaxonSet ts : allTs)
                    allTaxaOf.put(ts, ts.getTaxonSet().stream().map(Taxon::getID)
                        .collect(Collectors.toCollection(LinkedHashSet::new)));

                Map<TaxonSet, double[]> boundsOf = new HashMap<>();
                for (var ccp : ct.getCalibrationCladePriors())
                    boundsOf.put(ccp.getTaxa(), new double[]{ccp.getLower(), ccp.getUpper()});

                Map<TaxonSet, Integer> tsToIdx = new HashMap<>();
                for (int i = 0; i < allTs.size(); i++) tsToIdx.put(allTs.get(i), i);

                Map<TaxonSet, List<TaxonSet>> immChildren = new HashMap<>();
                for (TaxonSet parent : allTs) {
                    Set<String> pTaxa = allTaxaOf.get(parent);
                    List<TaxonSet> children = new ArrayList<>();
                    for (TaxonSet child : allTs) {
                        if (child == parent) continue;
                        Set<String> cTaxa = allTaxaOf.get(child);
                        if (!pTaxa.containsAll(cTaxa)) continue;
                        boolean hasIntermediate = allTs.stream()
                            .filter(mid -> mid != parent && mid != child)
                            .anyMatch(mid -> {
                                Set<String> mTaxa = allTaxaOf.get(mid);
                                return pTaxa.containsAll(mTaxa) && mTaxa.containsAll(cTaxa)
                                    && !mTaxa.equals(cTaxa) && !mTaxa.equals(pTaxa);
                            });
                        if (!hasIntermediate) children.add(child);
                    }
                    immChildren.put(parent, children);
                }

                // Build CalibrationEntry list preserving index ordering
                calibrationEntries.clear();
                colorIndex = 0;
                Map<TaxonSet, CalibrationEntry> tsToEntry = new HashMap<>();
                for (TaxonSet ts : allTs) {
                    String label = ts.getID() != null ? ts.getID()
                                 : "Clade_" + (tsToIdx.get(ts) + 1);
                    CalibrationEntry entry = new CalibrationEntry(label, COLORS[colorIndex % COLORS.length]);
                    colorIndex++;
                    double[] bounds = boundsOf.get(ts);
                    if (bounds != null) {
                        entry.lower = Double.isNaN(bounds[0]) ? null : bounds[0];
                        entry.upper = Double.isNaN(bounds[1]) ? null : bounds[1];
                    }
                    Set<String> childCoveredTaxa = immChildren.get(ts).stream()
                        .flatMap(child -> allTaxaOf.get(child).stream())
                        .collect(Collectors.toSet());
                    allTaxaOf.get(ts).stream()
                        .filter(t -> !childCoveredTaxa.contains(t) && taxa.contains(t))
                        .forEach(entry.directTaxa::add);
                    tsToEntry.put(ts, entry);
                    calibrationEntries.add(entry);
                }
                // Wire parent-child relationships
                for (TaxonSet ts : allTs) {
                    CalibrationEntry parent = tsToEntry.get(ts);
                    for (TaxonSet child : immChildren.get(ts))
                        parent.directChildCals.add(tsToEntry.get(child));
                }

                rebuildFn.run();
                dialog.close();
            } catch (Exception ex) {
                hint.setText("Parse error: " + ex.getMessage());
                hint.setStyle("-fx-font-size:11px; -fx-text-fill:#c00;");
            }
        });

        HBox btnRow = new HBox(8, okBtn, cancelBtn);
        btnRow.setPadding(new Insets(6, 0, 0, 0));
        VBox content = new VBox(8, ta, hint, btnRow);
        content.setPadding(new Insets(12));
        dialog.setScene(new Scene(content, 560, 220));
        dialog.showAndWait();
    }

    private void exportNewick() {
        Stage dialog = new Stage();
        dialog.initModality(Modality.APPLICATION_MODAL);
        dialog.setTitle("Export Constraint Tree (Newick)");
        TextArea ta = new TextArea(buildNewick());
        ta.setWrapText(true);
        ta.setPrefRowCount(6);
        ta.setEditable(false);
        Button closeBtn = new Button("Close");
        closeBtn.setOnAction(e -> dialog.close());
        VBox content = new VBox(8, ta, new HBox(closeBtn));
        content.setPadding(new Insets(12));
        dialog.setScene(new Scene(content, 560, 200));
        dialog.showAndWait();
    }

    private String buildNewick() {
        Set<CalibrationEntry> assignedAsChild = calibrationEntries.stream()
            .flatMap(e -> e.directChildCals.stream()).collect(Collectors.toSet());
        Set<String> ownedTaxa = calibrationEntries.stream()
            .flatMap(e -> e.directTaxa.stream()).collect(Collectors.toSet());

        List<String> parts = new ArrayList<>();
        for (CalibrationEntry e : calibrationEntries)
            if (!assignedAsChild.contains(e)) parts.add(buildNewickNode(e));
        for (String t : taxa)
            if (!ownedTaxa.contains(t)) parts.add(t);

        if (parts.isEmpty()) return "";
        if (parts.size() == 1) return parts.get(0) + ";";
        return "(" + String.join(",", parts) + ")[&virtualRoot=true];";
    }

    private String buildNewickNode(CalibrationEntry e) {
        List<String> kids = new ArrayList<>(e.directTaxa);
        for (CalibrationEntry child : e.directChildCals) kids.add(buildNewickNode(child));
        String ann = "[&name=" + e.label
            + (e.lower != null && e.upper != null ? ",lower=" + e.lower + ",upper=" + e.upper : "")
            + "]";
        if (kids.isEmpty()) return e.label;
        return "(" + String.join(",", kids) + ")" + ann;
    }

    // ── Save calibrations to model ────────────────────────────────────────────────

    private void saveCalibrations() {
        CalibratedBirthDeathSkylineModel cpp = (CalibratedBirthDeathSkylineModel) m_beastObject;
        String partition = partitionOf(cpp);

        cpp.calibrationsInput.get().clear();

        List<String> staleKeys = doc.pluginmap.keySet().stream()
            .filter(k -> (k.startsWith("CalibrationCladePrior.") || k.startsWith("TaxonSet."))
                      && k.endsWith("." + partition))
            .toList();
        staleKeys.forEach(doc.pluginmap::remove);

        for (CalibrationEntry entry : calibrationEntries) {
            List<String> allLeaves = allLeafTaxa(entry);
            if (allLeaves.isEmpty()) continue;

            List<Taxon> taxonList = new ArrayList<>();
            for (String t : allLeaves) {
                Taxon tx = new Taxon(t); tx.setID(t); taxonList.add(tx);
            }
            TaxonSet ts = new TaxonSet();
            ts.initByName("taxon", taxonList);
            ts.setID("TaxonSet." + entry.label + "." + partition);
            doc.addPlugin(ts);
            cpp.calibrationsInput.get().add(ts);

            if (entry.lower != null && entry.upper != null) {
                CalibrationCladePrior ccp = new CalibrationCladePrior();
                ccp.initByName("taxa", ts,
                    "lowerAge", new RealScalarParam<>(entry.lower, NonNegativeReal.INSTANCE),
                    "upperAge", new RealScalarParam<>(entry.upper, NonNegativeReal.INSTANCE));
                ccp.setID("CalibrationCladePrior." + entry.label + "." + partition);
                doc.addPlugin(ccp);
            }
        }
        sync();
    }

    private List<String> allLeafTaxa(CalibrationEntry e) {
        List<String> out = new ArrayList<>(e.directTaxa);
        for (CalibrationEntry child : e.directChildCals) {
            for (String t : allLeafTaxa(child))
                if (!out.contains(t)) out.add(t);
        }
        return out;
    }
}
