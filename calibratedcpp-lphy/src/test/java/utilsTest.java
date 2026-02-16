import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class utilsTest {
    @Test
    void testComputeParents() {
        Calibration calibration1 = new Calibration(new String[]{"1","2","3","4"});
        Calibration calibration2 = new Calibration(new String[]{"1","2"});
        Calibration calibration3 = new Calibration(new String[]{"3","4"});
        Calibration calibration4 = new Calibration(new String[]{"5","6"});
        Calibration calibration5 = new Calibration(new String[]{"root"});
        List<Calibration> calibrations = new ArrayList<>();
        calibrations.add(calibration5);
        calibrations.add(calibration1);
        calibrations.add(calibration2);
        calibrations.add(calibration3);
        calibrations.add(calibration4);

        boolean rootFlag=true;


        int[] observed = ConditionedMRCAPrior.computeParents(rootFlag, calibrations);
        int[] expected = new int[]{-1, 0,1,1,0};
        for (int i = 0; i < observed.length; i++) {
            assertEquals(expected[i], observed[i]);
            System.out.println(observed[i]);
            System.out.println(expected[i]);
        }
    }
}
