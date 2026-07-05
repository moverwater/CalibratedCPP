package calibratedcpp.beauti;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibration.ConstraintTree;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.stage.Modality;
import javafx.stage.Stage;

/**
 * Modal "Manage Calibrations" popup shared by all calibrated tree-prior editors. It edits a list of
 * {@link CalibrationEntry} objects in place (add/remove clades, assign taxa, nest child clades, set
 * bounds), offers Newick import/export of the constraint tree, and renders a live hierarchy preview
 * that updates on every edit. Persisting the entries back to the model is the caller's concern,
 * handled through the {@code onSaveClose} callback fired by "Save &amp; Close".
 */
class CalibrationManagerPanel {

    private final BeautiDoc doc;
    private final List<String> taxa;
    private final List<CalibrationEntry> entries;
    private final Supplier<String> nextColor;
    private final Runnable onSaveClose;

    private VBox previewBox; // live hierarchy preview; null while the popup is closed

    /**
     * @param doc         the BEAUti document (used only for Newick constraint-tree parsing)
     * @param taxa        the full taxon list of the alignment
     * @param entries     the mutable calibration list; edited in place so the caller sees changes
     * @param nextColor   supplies a colour for each newly-created clade
     * @param onSaveClose invoked when the user clicks "Save &amp; Close" (persist to the model here)
     */
    CalibrationManagerPanel(BeautiDoc doc, List<String> taxa, List<CalibrationEntry> entries,
                            Supplier<String> nextColor, Runnable onSaveClose) {
        this.doc = doc;
        this.taxa = taxa;
        this.entries = entries;
        this.nextColor = nextColor;
        this.onSaveClose = onSaveClose;
    }

    void show() {
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
            entries.add(new CalibrationEntry("Clade_" + (entries.size() + 1), nextColor.get()));
            rebuild.run();
        });
        importBtn.setOnAction(e -> importNewick(rebuild));
        exportBtn.setOnAction(e -> exportNewick());
        saveBtn.setOnAction(e -> { onSaveClose.run(); stage.close(); });

        // Preview pane: a nested outline of the calibration hierarchy that re-renders on every edit
        // (add/remove clade, rename, bounds, taxon or child assignment).
        Label previewTitle = new Label("Preview");
        previewTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 13px;");
        HBox previewHeader = new HBox(previewTitle);
        previewHeader.setPadding(new Insets(8, 10, 8, 10));
        previewHeader.setStyle("-fx-background-color: #ececec; -fx-border-color: #d0d0d0; -fx-border-width: 0 0 1 0;");
        previewBox = FXUtils.newVBox();
        previewBox.setSpacing(2);
        previewBox.setPadding(new Insets(10));
        ScrollPane previewScroll = new ScrollPane(previewBox);
        previewScroll.setFitToWidth(true);
        previewScroll.setStyle("-fx-background-color: white;");
        VBox previewPane = new VBox(previewHeader, previewScroll);
        previewPane.setPrefWidth(320);
        previewPane.setMinWidth(240);
        previewPane.setStyle("-fx-border-color: #d0d0d0; -fx-border-width: 0 0 0 1;");
        VBox.setVgrow(previewScroll, Priority.ALWAYS);

        rebuild.run();

        HBox middle = new HBox(scrollPane, previewPane);
        HBox.setHgrow(scrollPane, Priority.ALWAYS);
        VBox root = new VBox(topBar, middle, bottomBar);
        VBox.setVgrow(middle, Priority.ALWAYS);
        stage.setScene(new Scene(root, 1000, 620));
        stage.setOnHidden(e -> previewBox = null);
        stage.showAndWait();
    }

    // ── Card list ─────────────────────────────────────────────────────────────────

    private void rebuildCalibrationList(VBox cardBox, ScrollPane scrollPane) {
        double savedV = scrollPane.getVvalue();
        cardBox.getChildren().clear();
        for (CalibrationEntry entry : entries)
            cardBox.getChildren().add(buildCalibrationCard(entry, cardBox, scrollPane));
        refreshPreview();
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
        nameTf.textProperty().addListener((obs, o, n) -> {
            entry.label = n.replaceAll("[^\\w]", "_").isEmpty() ? entry.label : n.replaceAll("[^\\w]", "_");
            refreshPreview();
        });
        Button delBtn = new Button("×");
        delBtn.setStyle("-fx-background-color: #e74c3c; -fx-text-fill: white; -fx-padding: 2 7;");
        delBtn.setOnAction(e -> {
            for (CalibrationEntry other : entries) other.directChildCals.remove(entry);
            entries.remove(entry);
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
            refreshPreview();
        });
        boundsRow.getChildren().add(lowerTf);
        boundsRow.getChildren().add(new Label("Upper:"));
        TextField upperTf = boundsField(entry.upper);
        upperTf.textProperty().addListener((obs, o, n) -> {
            try { entry.upper = n.trim().isEmpty() ? null : Double.parseDouble(n); }
            catch (NumberFormatException ignored) {}
            refreshPreview();
        });
        boundsRow.getChildren().add(upperTf);

        card.getChildren().addAll(headerRow, boundsRow, new Separator());

        // Taxa section
        entry.taxaCbMap.clear();
        Label taxaHeading = new Label("Direct taxa:");
        taxaHeading.setStyle("-fx-font-size: 11px; -fx-font-weight: bold; -fx-text-fill: #555;");
        TextField taxaSearch = new TextField();
        taxaSearch.setPromptText("Search taxa…");
        taxaSearch.setStyle("-fx-font-size: 11px;");
        VBox taxaBox = new VBox(2);
        taxaBox.setPadding(new Insets(4, 4, 4, 8));
        for (String taxon : taxa) {
            boolean ownedByOther = entries.stream()
                .filter(e -> e != entry).anyMatch(e -> e.directTaxa.contains(taxon));
            CheckBox cb = new CheckBox(taxon);
            cb.setSelected(entry.directTaxa.contains(taxon));
            cb.setDisable(ownedByOther && !entry.directTaxa.contains(taxon));
            if (ownedByOther && !entry.directTaxa.contains(taxon))
                cb.setTooltip(new Tooltip("Assigned to " + entries.stream()
                    .filter(e -> e != entry && e.directTaxa.contains(taxon))
                    .map(e -> e.label).findFirst().orElse("?")));
            cb.selectedProperty().addListener((obs, o, n) -> {
                if (n) {
                    for (CalibrationEntry other : entries)
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
        // Live filter: show only checkboxes whose taxon label contains the query (case-insensitive).
        taxaSearch.textProperty().addListener((obs, o, n) -> {
            String q = n.trim().toLowerCase();
            for (Map.Entry<String, CheckBox> e : entry.taxaCbMap.entrySet()) {
                boolean visible = q.isEmpty() || e.getKey().toLowerCase().contains(q);
                e.getValue().setVisible(visible);
                e.getValue().setManaged(visible);
            }
        });
        card.getChildren().addAll(taxaHeading, taxaSearch, taxaSp);

        // Child calibrations section (only when there are others)
        entry.childCalCbMap.clear();
        List<CalibrationEntry> others = entries.stream().filter(e -> e != entry).toList();
        if (!others.isEmpty()) {
            Label calHeading = new Label("Child calibrations:");
            calHeading.setStyle("-fx-font-size: 11px; -fx-font-weight: bold; -fx-text-fill: #555;");
            VBox calBox = new VBox(2);
            calBox.setPadding(new Insets(4, 4, 4, 8));
            for (CalibrationEntry other : others) {
                boolean ownedByOther = entries.stream()
                    .filter(e -> e != entry).anyMatch(e -> e.directChildCals.contains(other));
                boolean wouldCycle = isAncestor(other, entry);
                CheckBox cb = new CheckBox(other.label);
                cb.setSelected(entry.directChildCals.contains(other));
                cb.setDisable((ownedByOther || wouldCycle) && !entry.directChildCals.contains(other));
                if (wouldCycle)
                    cb.setTooltip(new Tooltip("Would create a cycle"));
                else if (ownedByOther)
                    cb.setTooltip(new Tooltip("Child of " + entries.stream()
                        .filter(e -> e != entry && e.directChildCals.contains(other))
                        .map(e -> e.label).findFirst().orElse("?")));
                cb.selectedProperty().addListener((obs, o, n) -> {
                    if (n) {
                        for (CalibrationEntry e : entries)
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
        for (CalibrationEntry entry : entries) {
            for (String taxon : taxa) {
                CheckBox cb = entry.taxaCbMap.get(taxon);
                if (cb == null) continue;
                boolean owned = entry.directTaxa.contains(taxon);
                boolean ownedByOther = entries.stream()
                    .filter(e -> e != entry).anyMatch(e -> e.directTaxa.contains(taxon));
                cb.setDisable(ownedByOther && !owned);
                if (ownedByOther && !owned)
                    cb.setTooltip(new Tooltip("Assigned to " + entries.stream()
                        .filter(e -> e != entry && e.directTaxa.contains(taxon))
                        .map(e -> e.label).findFirst().orElse("?")));
                else
                    cb.setTooltip(null);
            }
            for (CalibrationEntry other : entries) {
                if (other == entry) continue;
                CheckBox cb = entry.childCalCbMap.get(other);
                if (cb == null) continue;
                boolean owned = entry.directChildCals.contains(other);
                boolean ownedByOther = entries.stream()
                    .filter(e -> e != entry).anyMatch(e -> e.directChildCals.contains(other));
                boolean wouldCycle = isAncestor(other, entry);
                cb.setDisable((ownedByOther || wouldCycle) && !owned);
                if (wouldCycle) cb.setTooltip(new Tooltip("Would create a cycle"));
                else if (ownedByOther)
                    cb.setTooltip(new Tooltip("Child of " + entries.stream()
                        .filter(e -> e != entry && e.directChildCals.contains(other))
                        .map(e -> e.label).findFirst().orElse("?")));
                else cb.setTooltip(null);
            }
        }
        refreshPreview();
    }

    private static TextField boundsField(Double value) {
        TextField tf = new TextField(value == null ? "" : String.valueOf(value));
        tf.setPrefWidth(90);
        tf.setPromptText("optional");
        return tf;
    }

    // ── Live preview ──────────────────────────────────────────────────────────────

    /** Re-renders the calibration hierarchy outline. Safe to call when the popup is closed. */
    private void refreshPreview() {
        if (previewBox == null) return;
        previewBox.getChildren().clear();

        Set<CalibrationEntry> assignedAsChild = entries.stream()
            .flatMap(e -> e.directChildCals.stream()).collect(Collectors.toSet());
        Set<String> ownedTaxa = entries.stream()
            .flatMap(e -> e.directTaxa.stream()).collect(Collectors.toSet());

        List<CalibrationEntry> roots = entries.stream()
            .filter(e -> !assignedAsChild.contains(e)).toList();
        List<String> freeTaxa = taxa.stream().filter(t -> !ownedTaxa.contains(t)).toList();

        if (roots.isEmpty() && freeTaxa.isEmpty()) {
            Label empty = new Label("No calibrations yet.");
            empty.setStyle("-fx-text-fill: #999; -fx-font-size: 11px;");
            previewBox.getChildren().add(empty);
            return;
        }

        for (CalibrationEntry root : roots)
            previewBox.getChildren().add(previewNode(root, 0));

        if (!freeTaxa.isEmpty()) {
            Label hdr = new Label("Unassigned taxa");
            hdr.setStyle("-fx-font-size: 10px; -fx-text-fill: #aaa;");
            hdr.setPadding(new Insets(6, 0, 0, 0));
            previewBox.getChildren().add(hdr);
            for (String t : freeTaxa)
                previewBox.getChildren().add(previewLeaf(t, 0, "#bbb"));
        }
    }

    private Node previewNode(CalibrationEntry e, int depth) {
        VBox box = new VBox(2);
        HBox header = new HBox(6);
        header.setAlignment(Pos.CENTER_LEFT);
        header.setPadding(new Insets(0, 0, 0, depth * 14.0));
        Region dot = new Region();
        dot.setMinSize(10, 10); dot.setMaxSize(10, 10);
        dot.setStyle("-fx-background-color: " + e.color + "; -fx-background-radius: 5;");
        Label name = new Label(e.label);
        name.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");
        header.getChildren().addAll(dot, name);
        if (e.lower != null && e.upper != null) {
            Label bounds = new Label("[" + fmtBound(e.lower) + ", " + fmtBound(e.upper) + "]");
            bounds.setStyle("-fx-font-size: 10px; -fx-text-fill: #2980b9;");
            header.getChildren().add(bounds);
        }
        box.getChildren().add(header);
        for (String t : e.directTaxa)
            box.getChildren().add(previewLeaf(t, depth + 1, "#666"));
        for (CalibrationEntry child : e.directChildCals)
            box.getChildren().add(previewNode(child, depth + 1));
        return box;
    }

    private static Node previewLeaf(String taxon, int depth, String color) {
        Label lbl = new Label("• " + taxon);
        lbl.setStyle("-fx-font-size: 11px; -fx-text-fill: " + color + ";");
        HBox row = new HBox(lbl);
        row.setPadding(new Insets(0, 0, 0, depth * 14.0 + 4));
        return row;
    }

    private static String fmtBound(double v) {
        return v == Math.floor(v) && !Double.isInfinite(v) ? String.valueOf((long) v) : String.valueOf(v);
    }

    // ── Newick import / export ──────────────────────────────────────────────────────

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
                entries.clear();
                Map<TaxonSet, CalibrationEntry> tsToEntry = new HashMap<>();
                for (TaxonSet ts : allTs) {
                    String label = ts.getID() != null ? ts.getID()
                                 : "Clade_" + (tsToIdx.get(ts) + 1);
                    CalibrationEntry entry = new CalibrationEntry(label, nextColor.get());
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
                    entries.add(entry);
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
        Set<CalibrationEntry> assignedAsChild = entries.stream()
            .flatMap(e -> e.directChildCals.stream()).collect(Collectors.toSet());
        Set<String> ownedTaxa = entries.stream()
            .flatMap(e -> e.directTaxa.stream()).collect(Collectors.toSet());

        List<String> parts = new ArrayList<>();
        for (CalibrationEntry e : entries)
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
}
