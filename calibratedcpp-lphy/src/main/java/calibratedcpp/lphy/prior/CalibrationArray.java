package calibratedcpp.lphy.prior;

import lphy.core.logger.TextFileFormatted;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

public class CalibrationArray implements TextFileFormatted {
    Calibration[] array;
    public CalibrationArray(Calibration[] array) {
        this.array = array;
    }

    public Calibration[] getCalibrationArray() {
        return array;
    }

    @Override
    public List<String> getTextForFile() {
        return List.of();
    }

    @Override
    public void writeToFile(BufferedWriter writer) {
        try {
            writer.write("calibrationTaxa\tage\n");

            for (int i = 0; i < array.length; i++) {
                StringBuilder sb = new StringBuilder();
                String[] taxa = array[i].getTaxa();
                double age = array[i].getAge();
                for (int j = 0; j < taxa.length; j++) {
                    if (j > 0) sb.append(",");
                    sb.append(taxa[j]);
                }

                sb.append("\t").append(age).append("\n");
                writer.write(sb.toString());
            }
            writer.flush();
        } catch (IOException e){
            e.printStackTrace();
        }
    }

    @Override
    public String getFileType() {
        return "_calibrations.log";
    }
}
