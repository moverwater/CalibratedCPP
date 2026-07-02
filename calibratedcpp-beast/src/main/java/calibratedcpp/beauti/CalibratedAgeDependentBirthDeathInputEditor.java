package calibratedcpp.beauti;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
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
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.Normal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.operator.RealRandomWalkOperator;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beastfx.app.beauti.TreeDistributionInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedAgeDependentBirthDeathModel;
import calibratedcpp.distribution.Erlang;
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
import javafx.scene.control.RadioButton;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.stage.Modality;
import javafx.stage.Stage;

/**
 * BEAUti panel for {@link CalibratedAgeDependentBirthDeathModel}. Shares the calibration-management
 * popup and root/origin conditioning row with {@link CalibratedCPPInputEditor} (duplicated rather
 * than factored out, to avoid touching that editor's delicate re-entrancy handling), but replaces
 * the skyline rate-parameterization UI with a much simpler birth rate / sampling probability /
 * lifetime-distribution picker, since this model has no time-varying rates.
 */
public class CalibratedAgeDependentBirthDeathInputEditor extends TreeDistributionInputEditor {

    // ── Calibration state ─────────────────────────────────────────────────────────

    static class CalibrationEntry {
        String label;
        final String color;
        final List<String> directTaxa = new ArrayList<>();
        final List<CalibrationEntry> directChildCals = new ArrayList<>();
        Double lower, upper;
        final Map<String, CheckBox> taxaCbMap = new HashMap<>();
        final Map<CalibrationEntry, CheckBox> childCalCbMap = new HashMap<>();
        CalibrationEntry(String label, String color) { this.label = label; this.color = color; }
    }

    private static final String[] COLORS = {
        "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
        "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac"
    };

    private List<String> taxa = List.of("Apple", "Banana", "Cherry");
    private final List<CalibrationEntry> calibrationEntries = new ArrayList<>();
    private int colorIndex = 0;

    private VBox lifetimeParamsBox;
    private RadioButton rootRb;
    private RadioButton originRb;

    private enum LifetimeKind {
        EXPONENTIAL("Exponential"), GAMMA("Gamma"), ERLANG("Erlang"), LOGNORMAL("Log-Normal");
        final String display;
        LifetimeKind(String display) { this.display = display; }
    }

    private enum ParamKind { POSITIVE, UNIT_INTERVAL, REAL }

    public CalibratedAgeDependentBirthDeathInputEditor(BeautiDoc doc) { super(doc); }
    public CalibratedAgeDependentBirthDeathInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibratedAgeDependentBirthDeathModel.class; }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int listItemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, listItemNr, isExpandOption, addButtons);

        CalibratedAgeDependentBirthDeathModel model = (CalibratedAgeDependentBirthDeathModel) m_beastObject;
        try { taxa = model.treeInput.get().getTaxonset().asStringList(); } catch (Exception ignored) {}

        if (calibrationEntries.isEmpty()) initCalibrationsFromModel(model);

        Button calibrationButton = new Button("Manage calibrations");
        calibrationButton.setOnAction(e -> openPopup());
        pane.getChildren().add(calibrationButton);

        String partition = partitionOf(model);

        // Birth rate
        if (doc.pluginmap.get("adBirthRate." + partition) instanceof RealScalarParam<?> birthRate) {
            bootstrapIfEstimated("adBirthRate", partition, birthRate, ParamKind.POSITIVE);
            pane.getChildren().add(scalarRow("Birth rate (λ)", "adBirthRate", partition, birthRate, ParamKind.POSITIVE));
        }

        // Lifetime distribution picker
        HBox distRow = FXUtils.newHBox();
        distRow.setSpacing(8);
        distRow.setPadding(new Insets(8, 0, 4, 0));
        distRow.setAlignment(Pos.CENTER_LEFT);
        distRow.getChildren().add(new Label("Lifetime distribution:"));
        ChoiceBox<String> distChoice = new ChoiceBox<>();
        for (LifetimeKind k : LifetimeKind.values()) distChoice.getItems().add(k.display);
        LifetimeKind currentKind = detectCurrentLifetimeKind(model);
        distChoice.getSelectionModel().select(currentKind.ordinal());
        distRow.getChildren().add(distChoice);
        pane.getChildren().add(distRow);

        // Bootstrap: ensure the currently-selected distribution's params/prior/operator exist
        // and that the model actually points at them (handles a freshly-applied template).
        Object activeDist = activateLifetimeDistribution(currentKind, partition);
        if (model.lifetimeDistributionInput.get() != activeDist) setLifetimeDistribution(model, activeDist);

        lifetimeParamsBox = FXUtils.newVBox();
        lifetimeParamsBox.setSpacing(10);
        pane.getChildren().add(lifetimeParamsBox);
        refreshLifetimeParamsBox(model);

        // Extant sampling probability (rho)
        if (doc.pluginmap.get("adRho." + partition) instanceof RealScalarParam<?> rho) {
            bootstrapIfEstimated("adRho", partition, rho, ParamKind.UNIT_INTERVAL);
            pane.getChildren().add(scalarRow("Extant sampling probability (ρ)", "adRho", partition, rho, ParamKind.UNIT_INTERVAL));
        }

        // Condition on root vs. set origin
        addConditionOnRootRow(pane, model);

        distChoice.getSelectionModel().selectedIndexProperty().addListener((obs, oldIdx, newIdx) -> {
            if (newIdx.intValue() < 0) return;
            LifetimeKind kind = LifetimeKind.values()[newIdx.intValue()];
            switchLifetimeDistribution(model, kind);
            refreshLifetimeParamsBox(model);
            sync();
        });
    }

    // ── Model initialization ──────────────────────────────────────────────────────

    private void initCalibrationsFromModel(CalibratedAgeDependentBirthDeathModel model) {
        String partition = partitionOf(model);
        List<TaxonSet> allTs = model.calibrationsInput.get();
        // Switching tree-prior templates builds a fresh model whose calibrations list is empty,
        // but BEAUti's removeSubNet only disconnects connectors — the TaxonSet.*/.partition
        // plugins survive in the pluginmap. Recover and re-attach them so calibrations persist.
        if (allTs.isEmpty()) recoverCalibrationTaxonSets(partition, allTs);
        if (allTs.isEmpty()) return;

        Map<TaxonSet, Set<String>> taxaOf = new LinkedHashMap<>();
        for (TaxonSet ts : allTs) {
            Set<String> tset = ts.getTaxonSet().stream().map(Taxon::getID)
                .collect(Collectors.toCollection(LinkedHashSet::new));
            taxaOf.put(ts, tset);
        }

        Map<TaxonSet, List<TaxonSet>> immChildren = new LinkedHashMap<>();
        for (TaxonSet parent : allTs) {
            Set<String> pTaxa = taxaOf.get(parent);
            List<TaxonSet> children = new ArrayList<>();
            for (TaxonSet child : allTs) {
                if (child == parent) continue;
                Set<String> cTaxa = taxaOf.get(child);
                if (!pTaxa.containsAll(cTaxa)) continue;
                boolean hasIntermediate = allTs.stream()
                    .filter(mid -> mid != parent && mid != child)
                    .anyMatch(mid -> {
                        Set<String> mTaxa = taxaOf.get(mid);
                        return pTaxa.containsAll(mTaxa) && mTaxa.containsAll(cTaxa)
                            && !mTaxa.equals(cTaxa) && !mTaxa.equals(pTaxa);
                    });
                if (!hasIntermediate) children.add(child);
            }
            immChildren.put(parent, children);
        }

        Map<TaxonSet, CalibrationEntry> tsToEntry = new LinkedHashMap<>();
        for (TaxonSet ts : allTs) {
            String label = tsLabelOf(ts.getID(), partition);
            CalibrationEntry entry = new CalibrationEntry(label, COLORS[colorIndex % COLORS.length]);
            colorIndex++;

            String ccpId = "CalibrationCladePrior." + label + "." + partition;
            BEASTInterface bi = doc.pluginmap.get(ccpId);
            if (bi instanceof CalibrationCladePrior ccp) {
                entry.lower = ccp.getLower();
                entry.upper = ccp.getUpper();
            }

            Set<String> childCoveredTaxa = immChildren.get(ts).stream()
                .flatMap(child -> taxaOf.get(child).stream())
                .collect(Collectors.toSet());
            taxaOf.get(ts).stream()
                .filter(t -> !childCoveredTaxa.contains(t) && taxa.contains(t))
                .forEach(entry.directTaxa::add);

            tsToEntry.put(ts, entry);
            calibrationEntries.add(entry);
        }

        for (TaxonSet ts : allTs) {
            CalibrationEntry parent = tsToEntry.get(ts);
            for (TaxonSet child : immChildren.get(ts))
                parent.directChildCals.add(tsToEntry.get(child));
        }
    }

    private static String tsLabelOf(String tsId, String partition) {
        if (tsId == null) return "unknown";
        String prefix = "TaxonSet.";
        String suffix = "." + partition;
        if (tsId.startsWith(prefix) && tsId.endsWith(suffix))
            return tsId.substring(prefix.length(), tsId.length() - suffix.length());
        return tsId;
    }

    /**
     * Re-attaches calibration TaxonSets that survive in the pluginmap (e.g. after a tree-prior
     * template switch) to the model's calibrations list. Only IDs of the form
     * {@code TaxonSet.{label}.{partition}} with a non-empty label are recovered, so the alignment's
     * own {@code TaxonSet.{partition}} and unrelated taxon sets are skipped.
     */
    private void recoverCalibrationTaxonSets(String partition, List<TaxonSet> target) {
        String prefix = "TaxonSet.";
        String suffix = "." + partition;
        int minLen = prefix.length() + suffix.length() + 1; // +1 ensures a non-empty label between
        for (BEASTInterface bi : doc.pluginmap.values()) {
            if (bi instanceof TaxonSet ts && ts.getID() != null
                    && ts.getID().startsWith(prefix) && ts.getID().endsWith(suffix)
                    && ts.getID().length() >= minLen
                    && !target.contains(ts)) {
                target.add(ts);
            }
        }
    }

    // ── Scalar parameter rows (birth rate / rho / lifetime-distribution params) ────

    private VBox scalarRow(String label, String name, String partition, RealScalarParam<?> scalar, ParamKind kind) {
        VBox box = FXUtils.newVBox();
        box.setSpacing(4);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");

        HBox headerRow = FXUtils.newHBox();
        headerRow.setSpacing(8);
        headerRow.setAlignment(Pos.CENTER_LEFT);
        headerRow.getChildren().add(new Label(label));

        boolean estimated = scalar instanceof StateNode sn && sn.isEstimatedInput.get();
        CheckBox estimateCb = new CheckBox("Estimate");
        estimateCb.setSelected(estimated);
        headerRow.getChildren().add(estimateCb);
        box.getChildren().add(headerRow);

        double initVal = scalar.valuesInput.get();
        TextField valueTf = compactField(initVal);
        HBox valRow = FXUtils.newHBox();
        valRow.setSpacing(4);
        valRow.getChildren().addAll(new Label("Value:"), valueTf);
        box.getChildren().add(valRow);

        valueTf.textProperty().addListener((obs, o, nw) -> {
            try {
                double v = Double.parseDouble(nw);
                @SuppressWarnings({"unchecked", "rawtypes"})
                RealScalarParam rsp = (RealScalarParam) scalar;
                rsp.valuesInput.setValue(v, rsp);
                rsp.initAndValidate();
            } catch (Exception ignored) {}
        });

        estimateCb.setOnAction(e -> {
            setEstimated(scalar, estimateCb.isSelected());
            if (estimateCb.isSelected()) ensurePriorAndOperator(kind, name, partition, scalar);
            sync();
        });

        return box;
    }

    private void bootstrapIfEstimated(String name, String partition, RealScalarParam<?> scalar, ParamKind kind) {
        if (scalar instanceof StateNode sn && sn.isEstimatedInput.get())
            ensurePriorAndOperator(kind, name, partition, scalar);
    }

    private void ensurePriorAndOperator(ParamKind kind, String name, String partition, RealScalarParam<?> scalar) {
        switch (kind) {
            case POSITIVE -> ensureLogNormalPriorAndScaler(name, partition, scalar);
            case UNIT_INTERVAL -> ensureUniformPriorAndScaler(name, partition, scalar);
            case REAL -> ensureNormalPriorAndRW(name, partition, scalar);
        }
    }

    private void ensureLogNormalPriorAndScaler(String name, String partition, RealScalarParam<?> scalar) {
        String priorId = name + ".prior." + partition;
        String scalerId = name + "Scaler." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                LogNormal prior = new LogNormal();
                prior.setInputValue("M", new RealScalarParam<>(0.0, Real.INSTANCE));
                prior.setInputValue("S", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));
                prior.setInputValue("param", scalar);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(scalerId)) {
            try {
                ScaleOperator scaler = new ScaleOperator();
                scaler.setInputValue("parameter", scalar);
                scaler.setInputValue("weight", 1.0);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    private void ensureUniformPriorAndScaler(String name, String partition, RealScalarParam<?> scalar) {
        String priorId = name + ".prior." + partition;
        String scalerId = name + "Scaler." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                Uniform prior = new Uniform();
                prior.setInputValue("lower", new RealScalarParam<>(0.0, Real.INSTANCE));
                prior.setInputValue("upper", new RealScalarParam<>(1.0, Real.INSTANCE));
                prior.setInputValue("param", scalar);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(scalerId)) {
            try {
                ScaleOperator scaler = new ScaleOperator();
                scaler.setInputValue("parameter", scalar);
                scaler.setInputValue("weight", 0.5);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    private void ensureNormalPriorAndRW(String name, String partition, RealScalarParam<?> scalar) {
        String priorId = name + ".prior." + partition;
        String rwId = name + "RW." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                Normal prior = new Normal();
                prior.setInputValue("mean", new RealScalarParam<>(0.0, Real.INSTANCE));
                prior.setInputValue("sigma", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));
                prior.setInputValue("param", scalar);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(rwId)) {
            try {
                RealRandomWalkOperator rw = new RealRandomWalkOperator();
                rw.setInputValue("scalar", scalar);
                rw.setInputValue("windowSize", 1.0);
                rw.setInputValue("weight", 1.0);
                rw.initAndValidate();
                pluginPut(rwId, rw);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    // ── Lifetime distribution ────────────────────────────────────────────────────

    private LifetimeKind detectCurrentLifetimeKind(CalibratedAgeDependentBirthDeathModel model) {
        Object d = model.lifetimeDistributionInput.get();
        if (d instanceof Erlang) return LifetimeKind.ERLANG;
        if (d instanceof Gamma) return LifetimeKind.GAMMA;
        if (d instanceof LogNormal) return LifetimeKind.LOGNORMAL;
        return LifetimeKind.EXPONENTIAL;
    }

    private Object activateLifetimeDistribution(LifetimeKind kind, String partition) {
        return switch (kind) {
            case EXPONENTIAL -> activateExponential(partition);
            case GAMMA -> activateGamma(partition);
            case ERLANG -> activateErlang(partition);
            case LOGNORMAL -> activateLogNormal(partition);
        };
    }

    private void switchLifetimeDistribution(CalibratedAgeDependentBirthDeathModel model, LifetimeKind kind) {
        String partition = partitionOf(model);
        deactivateCurrentLifetimeParams(model.lifetimeDistributionInput.get());
        setLifetimeDistribution(model, activateLifetimeDistribution(kind, partition));
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private void setLifetimeDistribution(CalibratedAgeDependentBirthDeathModel model, Object dist) {
        ((Input) model.lifetimeDistributionInput).setValue(dist, model);
    }

    private void deactivateCurrentLifetimeParams(Object dist) {
        if (dist instanceof Exponential e) {
            setEstimated(scalarOf(e.meanInput), false);
        } else if (dist instanceof Gamma g) {
            setEstimated(scalarOf(g.alphaInput), false);
            setEstimated(scalarOf(g.thetaInput), false);
        } else if (dist instanceof Erlang e) {
            setEstimated(scalarOf(e.scaleInput), false);
        } else if (dist instanceof LogNormal ln) {
            setEstimated(scalarOf(ln.MParameterInput), false);
            setEstimated(scalarOf(ln.SParameterInput), false);
        }
    }

    private Exponential activateExponential(String partition) {
        RealScalarParam<?> mean = ensureScalar("lifetimeMean", partition, 1.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = "lifetimeDistribution.exponential." + partition;
        if (doc.pluginmap.get(distId) instanceof Exponential existing) return existing;
        Exponential dist = new Exponential();
        dist.setInputValue("mean", mean);
        dist.initAndValidate();
        pluginPut(distId, dist);
        return dist;
    }

    private Gamma activateGamma(String partition) {
        RealScalarParam<?> shape = ensureScalar("lifetimeGammaShape", partition, 2.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        RealScalarParam<?> scale = ensureScalar("lifetimeGammaScale", partition, 1.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = "lifetimeDistribution.gamma." + partition;
        if (doc.pluginmap.get(distId) instanceof Gamma existing) return existing;
        Gamma dist = new Gamma();
        dist.setInputValue("alpha", shape);
        dist.setInputValue("theta", scale);
        dist.initAndValidate();
        pluginPut(distId, dist);
        return dist;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private Erlang activateErlang(String partition) {
        String shapeId = "lifetimeErlangShape." + partition;
        IntScalarParam shape;
        if (doc.pluginmap.get(shapeId) instanceof IntScalarParam existing) {
            shape = existing;
        } else {
            shape = new IntScalarParam(2, PositiveInt.INSTANCE);
            pluginPut(shapeId, shape);
        }
        RealScalarParam<?> scale = ensureScalar("lifetimeErlangScale", partition, 1.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = "lifetimeDistribution.erlang." + partition;
        if (doc.pluginmap.get(distId) instanceof Erlang existing) return existing;
        Erlang dist = new Erlang();
        dist.setInputValue("shape", shape);
        dist.setInputValue("scale", scale);
        dist.initAndValidate();
        pluginPut(distId, dist);
        return dist;
    }

    private LogNormal activateLogNormal(String partition) {
        RealScalarParam<?> m = ensureScalar("lifetimeLogNormalM", partition, 0.0, Real.INSTANCE, ParamKind.REAL);
        RealScalarParam<?> s = ensureScalar("lifetimeLogNormalS", partition, 1.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = "lifetimeDistribution.lognormal." + partition;
        if (doc.pluginmap.get(distId) instanceof LogNormal existing) return existing;
        LogNormal dist = new LogNormal();
        dist.setInputValue("M", m);
        dist.setInputValue("S", s);
        dist.initAndValidate();
        pluginPut(distId, dist);
        return dist;
    }

    /** Creates (or reuses) a scalar param, marking it estimated and ensuring its prior/operator exist. */
    private RealScalarParam<?> ensureScalar(String name, String partition, double defaultVal, Real domain, ParamKind kind) {
        String id = name + "." + partition;
        RealScalarParam<?> scalar;
        if (doc.pluginmap.get(id) instanceof RealScalarParam<?> existing) {
            scalar = existing;
        } else {
            RealScalarParam<Real> created = new RealScalarParam<>(defaultVal, domain);
            pluginPut(id, created);
            setEstimated(created, true);
            scalar = created;
        }
        if (scalar instanceof StateNode sn && sn.isEstimatedInput.get())
            ensurePriorAndOperator(kind, name, partition, scalar);
        return scalar;
    }

    private void refreshLifetimeParamsBox(CalibratedAgeDependentBirthDeathModel model) {
        lifetimeParamsBox.getChildren().clear();
        String partition = partitionOf(model);
        Object d = model.lifetimeDistributionInput.get();
        if (d instanceof Exponential e) {
            lifetimeParamsBox.getChildren().add(
                scalarRow("Mean lifetime", "lifetimeMean", partition, scalarOf(e.meanInput), ParamKind.POSITIVE));
        } else if (d instanceof Gamma g) {
            lifetimeParamsBox.getChildren().add(
                scalarRow("Shape (α)", "lifetimeGammaShape", partition, scalarOf(g.alphaInput), ParamKind.POSITIVE));
            lifetimeParamsBox.getChildren().add(
                scalarRow("Scale (θ)", "lifetimeGammaScale", partition, scalarOf(g.thetaInput), ParamKind.POSITIVE));
        } else if (d instanceof Erlang e) {
            lifetimeParamsBox.getChildren().add(erlangShapeRow(e));
            lifetimeParamsBox.getChildren().add(
                scalarRow("Scale (θ)", "lifetimeErlangScale", partition, scalarOf(e.scaleInput), ParamKind.POSITIVE));
        } else if (d instanceof LogNormal ln) {
            lifetimeParamsBox.getChildren().add(
                scalarRow("M (log mean)", "lifetimeLogNormalM", partition, scalarOf(ln.MParameterInput), ParamKind.REAL));
            lifetimeParamsBox.getChildren().add(
                scalarRow("S (log sd)", "lifetimeLogNormalS", partition, scalarOf(ln.SParameterInput), ParamKind.POSITIVE));
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private HBox erlangShapeRow(Erlang e) {
        HBox row = FXUtils.newHBox();
        row.setSpacing(6);
        row.setAlignment(Pos.CENTER_LEFT);
        row.getChildren().add(new Label("Shape (k, positive integer):"));
        IntScalarParam shape = (IntScalarParam) e.shapeInput.get();
        Spinner<Integer> spinner = new Spinner<>(1, 1000, (Integer) shape.valuesInput.get());
        spinner.setEditable(true);
        spinner.setPrefWidth(80);
        spinner.valueProperty().addListener((obs, o, n) -> {
            shape.valuesInput.setValue(n, shape);
            try { shape.initAndValidate(); } catch (Exception ignored) {}
        });
        row.getChildren().add(spinner);
        return row;
    }

    @SuppressWarnings("unchecked")
    private static RealScalarParam<?> scalarOf(Input<?> input) {
        return (input.get() instanceof RealScalarParam<?> r) ? r : null;
    }

    /** Registers a new plugin via doc.addPlugin (only call when the ID is not yet in pluginmap). */
    private void pluginPut(String id, BEASTInterface plugin) {
        plugin.setID(id);
        doc.addPlugin(plugin);
    }

    private static void setEstimated(RealScalarParam<?> scalar, boolean estimated) {
        if (scalar instanceof StateNode sn) sn.isEstimatedInput.setValue(estimated, sn);
    }

    private static String partitionOf(CalibratedAgeDependentBirthDeathModel model) {
        return model.getID().replaceFirst("^CalibratedAgeDependentBirthDeath\\.", "");
    }

    private static TextField compactField(double value) {
        TextField tf = new TextField(String.valueOf(value));
        tf.setPrefWidth(70);
        tf.setPadding(new Insets(2));
        return tf;
    }

    // ── Conditioning (root vs. origin) ──────────────────────────────────────────

    private void addConditionOnRootRow(Pane parent, CalibratedAgeDependentBirthDeathModel model) {
        String partition = partitionOf(model);

        VBox box = FXUtils.newVBox();
        box.setSpacing(6);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");

        Label heading = new Label("Conditioning");
        heading.setStyle("-fx-font-weight: bold;");
        box.getChildren().add(heading);

        ToggleGroup group = new ToggleGroup();
        rootRb   = new RadioButton("Condition on root");
        originRb = new RadioButton("Set origin time:");
        rootRb.setToggleGroup(group);
        originRb.setToggleGroup(group);

        boolean currentlyRoot = model.conditionOnRootInput.get();
        originRb.setDisable(hasFullTreeCalibration());
        rootRb.setSelected(currentlyRoot);
        originRb.setSelected(!currentlyRoot);

        String originParamId = "adOriginParam." + partition;
        RealScalarParam<?> originScalar =
            (doc.pluginmap.get(originParamId) instanceof RealScalarParam<?> rsp) ? rsp : null;

        if (currentlyRoot && model.originInput.get() != null)
            model.originInput.setValue(null, model);
        if (!currentlyRoot && originScalar != null && model.originInput.get() == null)
            model.originInput.setValue(originScalar, model);

        HBox originRow = FXUtils.newHBox();
        originRow.setSpacing(6);
        originRow.setAlignment(Pos.CENTER_LEFT);

        double initOriginVal = (originScalar != null) ? originScalar.valuesInput.get() : 1.0;
        TextField originTf = compactField(initOriginVal);
        originTf.setPrefWidth(90);

        boolean originEstimated = originScalar instanceof StateNode sn && sn.isEstimatedInput.get();
        CheckBox originEstimateCb = new CheckBox("Estimate");
        originEstimateCb.setSelected(originEstimated);

        originRow.getChildren().addAll(originTf, originEstimateCb);
        originRow.setVisible(!currentlyRoot);
        originRow.setManaged(!currentlyRoot);

        box.getChildren().addAll(rootRb, originRb, originRow);

        group.selectedToggleProperty().addListener((obs, oldTg, newTg) -> {
            if (newTg == null) return;
            boolean isRoot = (newTg == rootRb);
            model.conditionOnRootInput.setValue(isRoot, model);
            originRow.setVisible(!isRoot);
            originRow.setManaged(!isRoot);
            if (isRoot) {
                model.originInput.setValue(null, model);
            } else {
                RealScalarParam<?> os = (doc.pluginmap.get(originParamId) instanceof RealScalarParam<?> r) ? r : null;
                if (os != null) model.originInput.setValue(os, model);
            }
            sync();
        });

        if (originScalar != null) {
            originTf.textProperty().addListener((obs, o, nw) -> {
                try {
                    double v = Double.parseDouble(nw);
                    @SuppressWarnings({"unchecked", "rawtypes"})
                    RealScalarParam rsp = (RealScalarParam) originScalar;
                    rsp.valuesInput.setValue(v, rsp);
                    rsp.initAndValidate();
                } catch (Exception ignored) {}
            });
            originEstimateCb.setOnAction(e -> {
                setEstimated(originScalar, originEstimateCb.isSelected());
                if (originEstimateCb.isSelected())
                    ensurePriorAndOperator(ParamKind.POSITIVE, "adOriginParam", partition, originScalar);
                sync();
            });
        }

        parent.getChildren().add(box);
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
        CalibratedAgeDependentBirthDeathModel model = (CalibratedAgeDependentBirthDeathModel) m_beastObject;
        String partition = partitionOf(model);

        model.calibrationsInput.get().clear();

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
            model.calibrationsInput.get().add(ts);

            if (entry.lower != null && entry.upper != null) {
                CalibrationCladePrior ccp = new CalibrationCladePrior();
                ccp.initByName("taxa", ts,
                    "lowerAge", new RealScalarParam<>(entry.lower, beast.base.spec.domain.NonNegativeReal.INSTANCE),
                    "upperAge", new RealScalarParam<>(entry.upper, beast.base.spec.domain.NonNegativeReal.INSTANCE));
                ccp.setID("CalibrationCladePrior." + entry.label + "." + partition);
                doc.addPlugin(ccp);
            }
        }
        if (originRb != null) originRb.setDisable(hasFullTreeCalibration());
        sync();
    }

    private boolean hasFullTreeCalibration() {
        return calibrationEntries.stream()
            .anyMatch(e -> new HashSet<>(allLeafTaxa(e)).containsAll(taxa));
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
