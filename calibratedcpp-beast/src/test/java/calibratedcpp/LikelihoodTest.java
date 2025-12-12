package calibratedcpp;

import beast.base.evolution.speciation.CalibratedBirthDeathModel;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.tree.Tree;
import calibratedcpp.model.BirthDeathModel;
import calibration.CalibrationClade;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class LikelihoodTest {

    CalibratedCoalescentPointProcess cpp;
    BirthDeathModel birthDeathModel;
    CalibratedBirthDeathModel heled_and_drummond;
    Tree tree;
    List<CalibrationClade> calibrationsClades;
    List<CalibrationPoint> calibrationPoints;

    public LikelihoodTest() {
        cpp = new CalibratedCoalescentPointProcess();
        heled_and_drummond = new CalibratedBirthDeathModel();

        cpp.initByName("tree", tree,
                "treeModel", birthDeathModel,
                "calibrations", calibrationsClades);

        heled_and_drummond.initByName("tree", tree,
                "calibrations", calibrationPoints);
    }

    @Test
    public void calculateTreeLogLikelihoodTest() {
        assertEquals(cpp.calculateTreeLogLikelihood(tree),heled_and_drummond.calculateTreeLogLikelihood(tree),
                1e-6, "Calibrated CPP tree log likelihood does not math Heled and Drummond Tree log likelihood");
    }
}
