package calibratedcpp.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.Tensor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import calibratedcpp.SkylineParameter;
import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;

import java.util.ArrayList;
import java.util.List;

/**
 * BEAUti input editor for {@link SkylineParameter}.
 * Lets users set per-epoch rate values and (optionally) the change times
 * between epochs, together with the timesAreAges / timesAreRelative flags.
 */
public class SkylineParameterInputEditor extends InputEditor.Base {

    private SkylineParameter skylineParameter;
    private HBox changeTimesRow;
    private VBox changeTimesBox;
    private VBox valuesBox;
    private HBox valuesRow;

    public SkylineParameterInputEditor(BeautiDoc doc) {
        super(doc);
    }

    public SkylineParameterInputEditor() {
        super();
    }

    @Override
    public Class<?> type() {
        return SkylineParameter.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;
        pane = FXUtils.newVBox();

        skylineParameter = (SkylineParameter) input.get();

        addInputLabel();

        VBox mainBox = FXUtils.newVBox();
        mainBox.setSpacing(6);
        mainBox.setPadding(new Insets(4, 0, 4, 0));

        // ── Epoch count spinner ──────────────────────────────────────────────
        int nChanges = currentChangeCount();
        HBox epochRow = FXUtils.newHBox();
        epochRow.setSpacing(8);
        epochRow.getChildren().add(new Label("Number of epochs:"));
        Spinner<Integer> epochSpinner = new Spinner<>(1, 100, nChanges + 1);
        epochSpinner.setEditable(true);
        epochSpinner.setPrefWidth(70);
        epochRow.getChildren().add(epochSpinner);
        mainBox.getChildren().add(epochRow);

        // ── Change times area ────────────────────────────────────────────────
        changeTimesBox = FXUtils.newVBox();
        changeTimesBox.setSpacing(4);

        HBox timeFlagsRow = FXUtils.newHBox();
        timeFlagsRow.setSpacing(12);
        CheckBox timesAreAgesCb = new CheckBox("Times are ages");
        timesAreAgesCb.setTooltip(new Tooltip(
                "If checked, change times are ages relative to the most recent sample."));
        timesAreAgesCb.setSelected(skylineParameter.timesAreAgesInput.get());
        CheckBox timesAreRelativeCb = new CheckBox("Relative to process length (0–1)");
        timesAreRelativeCb.setSelected(skylineParameter.timesAreRelativeInput.get());
        timesAreRelativeCb.setTooltip(new Tooltip(
                "If checked, change times are expressed as fractions of the total process length (must be between 0 and 1)."));
        Button distributeEvenlyBtn = new Button("Distribute evenly");
        distributeEvenlyBtn.setTooltip(new Tooltip(
                "Space change times evenly over the process (fractions if relative, else 1, 2, ... )."));
        timeFlagsRow.getChildren().addAll(timesAreAgesCb, timesAreRelativeCb, distributeEvenlyBtn);
        changeTimesBox.getChildren().add(timeFlagsRow);

        changeTimesRow = FXUtils.newHBox();
        changeTimesRow.setSpacing(6);
        changeTimesBox.getChildren().add(changeTimesRow);
        mainBox.getChildren().add(changeTimesBox);

        // ── Values area ──────────────────────────────────────────────────────
        valuesBox = FXUtils.newVBox();
        valuesBox.setSpacing(4);
        valuesBox.getChildren().add(new Label("Epoch values (root → present):"));
        valuesRow = FXUtils.newHBox();
        valuesRow.setSpacing(6);
        valuesBox.getChildren().add(valuesRow);
        mainBox.getChildren().add(valuesBox);

        // initial render
        rebuildUI(nChanges, timesAreAgesCb, timesAreRelativeCb);

        // ── Listeners ────────────────────────────────────────────────────────
        epochSpinner.valueProperty().addListener((obs, oldVal, newVal) -> {
            int nNew = newVal - 1;
            ensureChangeTimes(nNew);
            ensureValues(newVal);
            rebuildUI(nNew, timesAreAgesCb, timesAreRelativeCb);
            sync();
        });

        timesAreAgesCb.selectedProperty().addListener((obs, o, n) -> {
            skylineParameter.timesAreAgesInput.setValue(n, skylineParameter);
            skylineParameter.initAndValidate();
            sync();
        });

        timesAreRelativeCb.selectedProperty().addListener((obs, o, n) -> {
            skylineParameter.timesAreRelativeInput.setValue(n, skylineParameter);
            // Redistribute evenly in the new scale so existing values don't violate constraints
            distributeEvenly(currentChangeCount(), n);
            rebuildUI(currentChangeCount(), timesAreAgesCb, timesAreRelativeCb);
            skylineParameter.initAndValidate();
            sync();
        });

        distributeEvenlyBtn.setOnAction(e -> {
            int nC = currentChangeCount();
            distributeEvenly(nC, timesAreRelativeCb.isSelected());
            rebuildUI(nC, timesAreAgesCb, timesAreRelativeCb);
            skylineParameter.initAndValidate();
            sync();
        });

        pane.getChildren().add(mainBox);
        getChildren().add(pane);
    }

    // ── UI rebuild ───────────────────────────────────────────────────────────

    private void rebuildUI(int nChanges,
                           CheckBox timesAreAgesCb, CheckBox timesAreRelativeCb) {
        boolean hasChanges = nChanges > 0;
        changeTimesBox.setVisible(hasChanges);
        changeTimesBox.setManaged(hasChanges);

        // Change-times text fields
        changeTimesRow.getChildren().clear();
        if (hasChanges) {
            RealVectorParam<?> ctParam = changeTimesParam();
            for (int i = 0; i < nChanges; i++) {
                double val = (ctParam != null && i < ctParam.size())
                        ? (double) ctParam.get(i) : 0.0;
                changeTimesRow.getChildren().add(epochLabel(i, nChanges,
                        timesAreAgesCb.isSelected()));
                TextField tf = timeField(val, i);
                changeTimesRow.getChildren().add(tf);
            }
        }

        // Values text fields
        valuesRow.getChildren().clear();
        int nEpochs = nChanges + 1;
        Tensor<?, Double> vParam = skylineParameter.valuesInput.get();
        for (int i = 0; i < nEpochs; i++) {
            double val = (vParam != null && i < vParam.size()) ? (double) vParam.get(i) : 1.0;
            valuesRow.getChildren().add(new Label("E" + (i + 1) + ":"));
            valuesRow.getChildren().add(valueField(val, i));
        }
    }

    private Label epochLabel(int i, int nChanges, boolean asAges) {
        String txt = asAges ? ("Age " + (i + 1) + ": ") : ("t" + (i + 1) + ": ");
        return new Label(txt);
    }

    // ── TextField factories ──────────────────────────────────────────────────

    private TextField timeField(double initial, int index) {
        TextField tf = compactField(initial);
        tf.textProperty().addListener((obs, o, n) -> {
            try {
                double v = Double.parseDouble(n);
                setChangeTime(index, v);
                skylineParameter.initAndValidate();
                sync();
            } catch (NumberFormatException ignored) {}
        });
        return tf;
    }

    private TextField valueField(double initial, int index) {
        TextField tf = compactField(initial);
        tf.textProperty().addListener((obs, o, n) -> {
            try {
                double v = Double.parseDouble(n);
                setEpochValue(index, v);
                skylineParameter.initAndValidate();
                sync();
            } catch (NumberFormatException ignored) {}
        });
        return tf;
    }

    private static TextField compactField(double value) {
        TextField tf = new TextField(String.valueOf(value));
        tf.setPrefWidth(65);
        tf.setPadding(new Insets(2));
        HBox.setMargin(tf, new Insets(0, 4, 0, 0));
        return tf;
    }

    // ── Parameter mutation ───────────────────────────────────────────────────

    @SuppressWarnings("unchecked")
    private void ensureChangeTimes(int nChanges) {
        if (nChanges == 0) {
            skylineParameter.changeTimesInput.setValue(null, skylineParameter);
            return;
        }
        boolean relative = skylineParameter.timesAreRelativeInput.get();
        RealVectorParam<NonNegativeReal> ctParam = changeTimesParam();
        if (ctParam == null) {
            double[] init = evenlySpaced(nChanges, relative);
            ctParam = new RealVectorParam<>(init, NonNegativeReal.INSTANCE);
            ctParam.setID(skylineParameter.getID() + ".changeTimes");
            doc.addPlugin(ctParam);
            skylineParameter.changeTimesInput.setValue(ctParam, skylineParameter);
        } else {
            double[] old = toDoubles(ctParam);
            double[] updated = new double[nChanges];
            for (int i = 0; i < nChanges; i++)
                updated[i] = i < old.length ? old[i] : evenlySpaced(nChanges, relative)[i];
            resetVectorValues(ctParam, updated);
        }
    }

    private void distributeEvenly(int nChanges, boolean relative) {
        RealVectorParam<?> p = changeTimesParam();
        if (p == null || nChanges == 0) return;
        resetVectorValues(p, evenlySpaced(nChanges, relative));
    }

    private static double[] evenlySpaced(int n, boolean relative) {
        double[] v = new double[n];
        for (int i = 0; i < n; i++)
            v[i] = relative ? (double)(i + 1) / (n + 1) : (double)(i + 1);
        return v;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private void ensureValues(int nEpochs) {
        Tensor<?, Double> existing = skylineParameter.valuesInput.get();
        Real domain = inferDomain(existing);

        if (existing instanceof RealVectorParam<?> vp) {
            double[] old = toDoubles(vp);
            double[] updated = new double[nEpochs];
            for (int i = 0; i < nEpochs; i++)
                updated[i] = i < old.length ? old[i] : (old.length > 0 ? old[old.length - 1] : 1.0);
            resetVectorValues(vp, updated);
        } else {
            // was a scalar — promote to vector
            double scalar = (existing != null) ? (double) existing.get(0) : 1.0;
            double[] init = new double[nEpochs];
            for (int i = 0; i < nEpochs; i++) init[i] = scalar;
            RealVectorParam<?> vp = new RealVectorParam(init, domain);
            vp.setID(skylineParameter.getID() + ".values");
            doc.addPlugin(vp);
            skylineParameter.valuesInput.setValue(vp, skylineParameter);
        }
    }

    private void setChangeTime(int index, double value) {
        RealVectorParam<?> p = changeTimesParam();
        if (p == null) return;
        double[] vals = toDoubles(p);
        if (index < vals.length) {
            vals[index] = value;
            resetVectorValues(p, vals);
        }
    }

    private void setEpochValue(int index, double value) {
        Tensor<?, Double> t = skylineParameter.valuesInput.get();
        if (!(t instanceof RealVectorParam<?> vp)) return;
        double[] vals = toDoubles(vp);
        if (index < vals.length) {
            vals[index] = value;
            resetVectorValues(vp, vals);
        }
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private int currentChangeCount() {
        Tensor<?, ?> ct = skylineParameter.changeTimesInput.get();
        return ct == null ? 0 : ct.size();
    }

    @SuppressWarnings("unchecked")
    private RealVectorParam<NonNegativeReal> changeTimesParam() {
        Tensor<?, ?> t = skylineParameter.changeTimesInput.get();
        return (t instanceof RealVectorParam<?>) ? (RealVectorParam<NonNegativeReal>) t : null;
    }

    private static double[] toDoubles(RealVectorParam<?> vp) {
        double[] out = new double[vp.size()];
        for (int i = 0; i < out.length; i++) out[i] = (double) vp.get(i);
        return out;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static void resetVectorValues(RealVectorParam<?> vp, double[] values) {
        // BEAST2's Input.setValue(List, ...) appends rather than replaces — mutate in-place.
        List<Double> existing = ((RealVectorParam<Real>) vp).valuesInput.get();
        existing.clear();
        for (double v : values) existing.add(v);
        // initAndValidate() uses max(dimensionInput, list.size()); reset dimensionInput first
        // so a previously-larger cached size doesn't re-expand the array after a shrink.
        vp.dimensionInput.setValue(values.length, vp);
        vp.initAndValidate();
    }

    private static Real inferDomain(Tensor<?, Double> t) {
        if (t instanceof RealScalarParam<?> sp) return sp.domainTypeInput.get();
        if (t instanceof RealVectorParam<?> vp) return vp.domainTypeInput.get();
        return Real.INSTANCE;
    }
}
