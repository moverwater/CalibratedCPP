import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import static calibratedcpp.lphy.prior.ConditionedMRCAPrior.mapBetaNodes;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class priorTest {
    @Test
    void testMapBetaArray_nonOverlap() {
        Double[] upperBounds = new Double[]{7.0,4.5, 3.0, 3.0, 4.0};
        Double[] lowerBounds = new Double[]{5.0, 4.0, 2.0, 1.0, 2.0};
        int[] parent = new int[]{-1,0,1,1,0};
        boolean rootFlag = true;

        boolean[] expect = new boolean[]{false,false,false,false,false};
        boolean[] actual = mapBetaNodes(parent.length, rootFlag,parent,upperBounds,lowerBounds);

        for (int i = 0; i < actual.length; i++) {
//            System.out.println(actual[i]);
//            System.out.println(expect[i]);
            assertEquals(expect[i],actual[i]);
        }
    }

    @Test
    void testMapBetaArray_withOverlap() {
        Double[] upperBounds = new Double[]{7.0,4.5, 3.0, 4.1, 4.0};
        Double[] lowerBounds = new Double[]{5.0, 4.0, 2.0, 1.0, 2.0};
        int[] parent = new int[]{-1,0,1,1,0};
        boolean rootFlag = true;

        boolean[] expect = new boolean[]{false,false,false,true,false};
        boolean[] actual = mapBetaNodes(parent.length, rootFlag,parent,upperBounds,lowerBounds);

        for (int i = 0; i < actual.length; i++) {
//            System.out.println("actual: " + actual[i]);
//            System.out.println("expect: " + expect[i]);
            assertEquals(expect[i],actual[i]);
        }
    }

    @Test
    void outputArray() {
        String[][] calibrationNames = new String[][]{new String[]{"1","2","3","4","5"}, new String[]{"1","2","3","4"}, new String[]{"1","2"}, new String[]{"3","4"}};
        Double[] upperBounds = new Double[]{6.0, 4.0, 3.0,2.5};
        Double[] lowerBounds = new Double[]{5.0,2.4, 2.0, 1.0};
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(new Value<>("", calibrationNames), new Value<>("", true),
                new Value<>("", upperBounds), new Value<>("", lowerBounds), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();
        assertEquals(observed.length,calibrationNames.length);
        for (int i = 0; i < observed.length; i++) {
            if (i == 0){
                assertEquals(5, observed[i].getTaxa().length);
                assertTrue (observed[i].getAge() > 5 && observed[i].getAge() < 6);
            } else {
                assertEquals(observed[i].getTaxa().length,calibrationNames[i].length);
                if (i == 1){
                    assertTrue (observed[i].getAge() > 2.4 && observed[i].getAge() < 4);
                } else if (i == 2){
                    assertTrue (observed[i].getAge() > 2 && observed[i].getAge() < observed[i-1].getAge());
                } else {
                    assertTrue (observed[i].getAge() > 1 && observed[i].getAge() < 2.5);
                }
            }
        }
    }

    @Test
    void outputArrayNoRoot() {
        String[][] calibrationNames = new String[][]{new String[]{"1","2","3","4"}, new String[]{"1","2"}, new String[]{"3","4"}};
        Double[] upperBounds = new Double[]{4.0, 3.0,2.5};
        Double[] lowerBounds = new Double[]{2.4, 2.0, 1.0};
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(new Value<>("", calibrationNames), new Value<>("", false),
                new Value<>("", upperBounds), new Value<>("", lowerBounds), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();
        assertEquals(observed.length,calibrationNames.length);
        for (int i = 0; i < observed.length; i++) {
            assertEquals(observed[i].getTaxa().length,calibrationNames[i].length);
            if (i == 0){
                assertTrue (observed[i].getAge() > 2.4 && observed[i].getAge() < 4);
            } else if (i == 1){
                assertTrue (observed[i].getAge() > 2 && observed[i].getAge() < observed[i-1].getAge());
            } else {
                assertTrue (observed[i].getAge() > 1 && observed[i].getAge() < 2.5);
            }
        }
    }

    @Test
    void inferenceTest() {
        String[][] calibrationNames = new String[][]{new String[]{"1","2","3"}, new String[]{"1","2","3","4"}, new String[]{"6","7"}};
        Double[] upperBounds = new Double[]{1.5, 1.9, 1.5};
        Double[] lowerBounds = new Double[]{1.0, 1.5, 1.1};
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(new Value<>("", calibrationNames), new Value<>("", false),
                new Value<>("", upperBounds), new Value<>("", lowerBounds), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();

        for (int i = 0; i < observed.length; i++) {
            assertEquals(observed[i].getTaxa().length,calibrationNames[i].length);
            if (i == 0){
                System.out.println(observed[i].getAge());
                assertTrue (observed[i].getAge() >1 && observed[i].getAge() < 1.5);
            } else if (i == 1){
                assertTrue (observed[i].getAge() < 1.9 && observed[i].getAge() > observed[i-1].getAge());
            } else {
                assertTrue (observed[i].getAge() > 1.1 && observed[i].getAge() < 1.5);
            }
        }
    }
}
