package calibratedcpp.beauti;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javafx.scene.control.CheckBox;

/**
 * In-memory model of one calibration clade edited in the {@link CalibrationManagerPanel}.
 * A clade owns a set of direct taxa and may nest child clades, so a collection of entries
 * forms a forest. The {@code *CbMap} fields hold transient references to the checkboxes that
 * represent this entry in the currently-open card, used for in-place enable/disable updates.
 */
class CalibrationEntry {
    String label;
    final String color;
    final List<String> directTaxa = new ArrayList<>();
    final List<CalibrationEntry> directChildCals = new ArrayList<>();
    Double lower, upper;
    final Map<String, CheckBox> taxaCbMap = new HashMap<>();
    final Map<CalibrationEntry, CheckBox> childCalCbMap = new HashMap<>();

    CalibrationEntry(String label, String color) {
        this.label = label;
        this.color = color;
    }

    /** All leaf taxa under this clade, including those contributed by nested child clades. */
    List<String> allLeafTaxa() {
        List<String> out = new ArrayList<>(directTaxa);
        for (CalibrationEntry child : directChildCals)
            for (String t : child.allLeafTaxa())
                if (!out.contains(t)) out.add(t);
        return out;
    }
}
