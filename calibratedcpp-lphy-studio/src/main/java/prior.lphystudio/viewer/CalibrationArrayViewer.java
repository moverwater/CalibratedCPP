package prior.lphystudio.viewer;

import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;

import javax.swing.*;

public class CalibrationArrayViewer implements Viewer {
    public CalibrationArrayViewer() {}

    @Override
    public boolean match(Object value) {
        return value instanceof CalibrationArray||
                (value instanceof Value && ((Value) value).value() instanceof CalibrationArray);
    }

    @Override
    public JComponent getViewer(Object value) {
        if (match(value)) {
            return new JTextArea(value.toString());
        }
        String text = ((Value<CalibrationArray>) value).value().toString();
        return new JTextArea(text);
    }

    @Override
    public String toString() {
        return "Calibration Array Viewer";
    }
}
