package calibratedcpp.beauti;

import beast.base.core.Input;
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
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedAgeDependentBirthDeathModel;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.distribution.Erlang;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;

/**
 * BEAUti panel for {@link CalibratedAgeDependentBirthDeathModel}. Inherits the
 * calibration-management popup, live preview, conditioning row, and model persistence from
 * {@link CalibratedCPPInputEditor}; this class contributes only the model-specific
 * parameter UI: birth rate, extant sampling probability, and a lifetime-distribution picker
 * (exponential / gamma / Erlang / log-normal) with its associated parameter editors.
 */
public class CalibratedAgeDependentBirthDeathInputEditor extends CalibratedCPPInputEditor {

    private enum LifetimeKind {
        EXPONENTIAL("Exponential"), GAMMA("Gamma"), ERLANG("Erlang"), LOGNORMAL("Log-Normal");
        final String display;
        LifetimeKind(String display) { this.display = display; }
    }

    private enum ParamKind { POSITIVE, UNIT_INTERVAL, REAL }

    private VBox lifetimeParamsBox;

    public CalibratedAgeDependentBirthDeathInputEditor(BeautiDoc doc) { super(doc); }
    public CalibratedAgeDependentBirthDeathInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibratedAgeDependentBirthDeathModel.class; }

    // ── Base-class hooks ────────────────────────────────────────────────────────────

    @Override
    protected String modelIdPrefix() { return "CalibratedAgeDependentBirthDeath"; }

    @Override
    protected String originParamId(String partition) { return "adOriginParam." + partition; }

    @Override
    protected void ensureOriginPriorAndOperator(String partition, RealScalarParam<?> scalar) {
        ensurePriorAndOperator(ParamKind.POSITIVE, "adOriginParam", partition, scalar);
    }

    @Override
    protected void buildModelUI(Pane pane, CalibratedCoalescentPointProcess baseModel) {
        CalibratedAgeDependentBirthDeathModel model = (CalibratedAgeDependentBirthDeathModel) baseModel;
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

        distChoice.getSelectionModel().selectedIndexProperty().addListener((obs, oldIdx, newIdx) -> {
            if (newIdx.intValue() < 0) return;
            LifetimeKind kind = LifetimeKind.values()[newIdx.intValue()];
            switchLifetimeDistribution(model, kind);
            refreshLifetimeParamsBox(model);
            sync();
        });
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
}
