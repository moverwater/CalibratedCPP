package calibratedcpp.beauti;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beastfx.app.beauti.TreeDistributionInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibrationprior.CalibrationCladePrior;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;

/**
 * Shared BEAUti editor scaffolding for calibrated tree priors that extend
 * {@link CalibratedCoalescentPointProcess}. It owns the calibration hierarchy state, the
 * "Manage calibrations" button (which opens a shared {@link CalibrationManagerPanel}), the
 * conditioning (root vs. origin) row, and persistence of calibrations to the model. Concrete
 * subclasses supply only the model-specific parameter UI and a handful of naming/estimation hooks.
 */
public abstract class CalibratedCPPInputEditor extends TreeDistributionInputEditor {

    protected static final String[] COLORS = {
        "#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f",
        "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac"
    };

    protected List<String> taxa = List.of("Apple", "Banana", "Cherry");
    protected final List<CalibrationEntry> calibrationEntries = new ArrayList<>();
    protected int colorIndex = 0;

    protected RadioButton rootRb;
    protected RadioButton originRb;

    protected CalibratedCPPInputEditor(BeautiDoc doc) { super(doc); }
    protected CalibratedCPPInputEditor() { super(); }

    // ── Model-specific hooks ───────────────────────────────────────────────────────

    /** The model's ID prefix, e.g. {@code "CalibratedCPP"} — used to derive the partition name. */
    protected abstract String modelIdPrefix();

    /** pluginmap ID of this partition's origin-time scalar, e.g. {@code "originParam." + partition}. */
    protected abstract String originParamId(String partition);

    /** Builds the model-specific parameter UI (rate parameterization, lifetime distribution, …). */
    protected abstract void buildModelUI(Pane pane, CalibratedCoalescentPointProcess model);

    /** Ensures a prior and operator exist for the origin scalar once its Estimate box is ticked. */
    protected abstract void ensureOriginPriorAndOperator(String partition, RealScalarParam<?> scalar);

    /**
     * Hook invoked after the origin has just been reconnected (switching to origin-time mode).
     * Default is a no-op; the skyline editor overrides it to mark rate change-times as relative.
     */
    protected void onOriginReconnected(CalibratedCoalescentPointProcess model) {}

    // ── Shared init template ────────────────────────────────────────────────────────

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int listItemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, listItemNr, isExpandOption, addButtons);

        CalibratedCoalescentPointProcess model = (CalibratedCoalescentPointProcess) m_beastObject;
        try { taxa = model.treeInput.get().getTaxonset().asStringList(); } catch (Exception ignored) {}

        if (calibrationEntries.isEmpty()) initCalibrationsFromModel(model);

        Button calibrationButton = new Button("Manage calibrations");
        calibrationButton.setOnAction(e -> new CalibrationManagerPanel(
            doc, taxa, calibrationEntries,
            () -> COLORS[colorIndex++ % COLORS.length],
            this::saveCalibrations).show());
        pane.getChildren().add(calibrationButton);

        buildModelUI(pane, model);
        addConditionOnRootRow(pane, model);
    }

    // ── Model initialization ──────────────────────────────────────────────────────

    /** Rebuilds the in-memory calibration forest from the model's TaxonSets and their CCP bounds. */
    protected void initCalibrationsFromModel(CalibratedCoalescentPointProcess model) {
        String partition = partitionOf(model);
        // calibrationsInput is the authoritative list — it includes all TaxonSets, even entries with
        // no bounds (and therefore no CCP in the pluginmap).
        List<TaxonSet> allTs = model.calibrationsInput.get();
        // Switching tree-prior templates builds a fresh model whose calibrations list is empty, but
        // BEAUti's removeSubNet only disconnects connectors — the TaxonSet.*/.partition plugins
        // survive in the pluginmap. Recover and re-attach them so calibrations persist.
        if (allTs.isEmpty()) recoverCalibrationTaxonSets(partition, allTs);
        if (allTs.isEmpty()) return;

        // Precompute leaf-taxa sets
        Map<TaxonSet, Set<String>> taxaOf = new LinkedHashMap<>();
        for (TaxonSet ts : allTs) {
            Set<String> tset = ts.getTaxonSet().stream().map(Taxon::getID)
                .collect(Collectors.toCollection(LinkedHashSet::new));
            taxaOf.put(ts, tset);
        }

        // Build immediate-children map (same containment logic as Newick import)
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

        // Create CalibrationEntry objects
        Map<TaxonSet, CalibrationEntry> tsToEntry = new LinkedHashMap<>();
        for (TaxonSet ts : allTs) {
            String label = tsLabelOf(ts.getID(), partition);
            CalibrationEntry entry = new CalibrationEntry(label, COLORS[colorIndex % COLORS.length]);
            colorIndex++;

            // Look up bounds from the corresponding CCP (may not exist if no bounds were set)
            String ccpId = "CalibrationCladePrior." + label + "." + partition;
            BEASTInterface bi = doc.pluginmap.get(ccpId);
            if (bi instanceof CalibrationCladePrior ccp) {
                entry.lower = ccp.getLower();
                entry.upper = ccp.getUpper();
            }

            // directTaxa = this entry's leaf taxa minus those covered by immediate children
            Set<String> childCoveredTaxa = immChildren.get(ts).stream()
                .flatMap(child -> taxaOf.get(child).stream())
                .collect(Collectors.toSet());
            taxaOf.get(ts).stream()
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
    }

    protected static String tsLabelOf(String tsId, String partition) {
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
    protected void recoverCalibrationTaxonSets(String partition, List<TaxonSet> target) {
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

    // ── Conditioning (root vs. origin) row ──────────────────────────────────────────

    protected void addConditionOnRootRow(Pane parent, CalibratedCoalescentPointProcess model) {
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

        String originParamId = originParamId(partition);
        RealScalarParam<?> originScalar =
            (doc.pluginmap.get(originParamId) instanceof RealScalarParam<?> rsp) ? rsp : null;

        // Fix up stale state: old XML may have origin connected even when conditionOnRoot=true
        // (e.g. saved before this fix, or from a template that always wired origin).
        if (currentlyRoot && model.originInput.get() != null)
            model.originInput.setValue(null, model);
        // Symmetrically, if we're in origin mode but origin is not connected, connect it now.
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
                // Disconnect origin from the model so it is absent from the XML.
                model.originInput.setValue(null, model);
            } else {
                // Reconnect origin, then let the subclass react (e.g. mark change times relative).
                RealScalarParam<?> os = (doc.pluginmap.get(originParamId) instanceof RealScalarParam<?> r) ? r : null;
                if (os != null) model.originInput.setValue(os, model);
                onOriginReconnected(model);
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
                if (originEstimateCb.isSelected()) ensureOriginPriorAndOperator(partition, originScalar);
                sync();
            });
        }

        parent.getChildren().add(box);
    }

    // ── Save calibrations to model ──────────────────────────────────────────────────

    protected void saveCalibrations() {
        CalibratedCoalescentPointProcess model = (CalibratedCoalescentPointProcess) m_beastObject;
        String partition = partitionOf(model);

        model.calibrationsInput.get().clear();

        List<String> staleKeys = doc.pluginmap.keySet().stream()
            .filter(k -> (k.startsWith("CalibrationCladePrior.") || k.startsWith("TaxonSet."))
                      && k.endsWith("." + partition))
            .toList();
        staleKeys.forEach(doc.pluginmap::remove);

        for (CalibrationEntry entry : calibrationEntries) {
            List<String> allLeaves = entry.allLeafTaxa();
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
                    "lowerAge", new RealScalarParam<>(entry.lower, NonNegativeReal.INSTANCE),
                    "upperAge", new RealScalarParam<>(entry.upper, NonNegativeReal.INSTANCE));
                ccp.setID("CalibrationCladePrior." + entry.label + "." + partition);
                doc.addPlugin(ccp);
            }
        }
        if (originRb != null) originRb.setDisable(hasFullTreeCalibration());
        sync();
    }

    protected boolean hasFullTreeCalibration() {
        return calibrationEntries.stream()
            .anyMatch(e -> new HashSet<>(e.allLeafTaxa()).containsAll(taxa));
    }

    // ── Shared helpers ──────────────────────────────────────────────────────────────

    protected String partitionOf(CalibratedCoalescentPointProcess model) {
        return model.getID().replaceFirst("^" + Pattern.quote(modelIdPrefix()) + "\\.", "");
    }

    protected void pluginPut(String id, BEASTInterface plugin) {
        plugin.setID(id);
        doc.addPlugin(plugin);
    }

    protected static void setEstimated(Object node, boolean estimated) {
        if (node instanceof StateNode sn) sn.isEstimatedInput.setValue(estimated, sn);
    }

    protected static TextField compactField(double value) {
        TextField tf = new TextField(String.valueOf(value));
        tf.setPrefWidth(70);
        tf.setPadding(new Insets(2));
        return tf;
    }
}
