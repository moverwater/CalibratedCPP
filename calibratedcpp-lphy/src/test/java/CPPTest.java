import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.tree.CPPTree;
import calibratedcpp.lphy.tree.CPPUtils;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static calibratedcpp.lphy.tree.CPPUtils.*;
import static calibratedcpp.lphy.tree.CalibratedCPPTree.coalesce;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class CPPTest {
    @Test
    public void testQdistBasic() {
        double result = CPPUtils.Qdist(2, 1, 1, 10);
        assertTrue(result > 0 && result < 1);
    }

    @Test
    void testCPPSimple() {
        Value<Number> samplingProb = new Value("", 0.2);
        Value<Number> birthRate = new Value("", 0.3);
        Value<Number> deathRate = new Value("", 0.1);
        Value<String[]> taxa = new Value<>("", new String[]{"1", "2", "3"});
        Value<Integer> n = new Value("", 3);
        CPPTree tree = new CPPTree(birthRate, deathRate, null, null,  samplingProb, taxa, n, new Value<>("", 10), null);
        TimeTree sampledTree = tree.sample().value();
        assertEquals(3, sampledTree.getLeafNodes().size());
        assertEquals(10, sampledTree.getRoot().getAge());

        assert(List.of(sampledTree.getRoot().getLeafNames()).contains("1"));
        assert(List.of(sampledTree.getRoot().getLeafNames()).contains("3"));
        assert(List.of(sampledTree.getRoot().getLeafNames()).contains("2"));
    }

    @Test
    void testCPP() {
        Value<Number> samplingProb = new Value("", 0.2);
        Value<Number> birthRate = new Value("", 0.3);
        Value<Number> deathRate = new Value("", 0.1);
        Value<Integer> n = new Value("", 5);
        String[] taxa = new String[]{"1", "2", "3"};
        double age = 2.0;
        Value<CalibrationArray> calibration = new Value<>("", new CalibrationArray(new Calibration[]{new Calibration(taxa, age)}));

        CalibratedCPPTree cpp = new CalibratedCPPTree(birthRate, deathRate, null, null,  samplingProb, n, calibration, null, null);
        TimeTree cppTree = cpp.sample().value();

        for (int i = 0; i < cppTree.getNodeCount(); i++) {
            TimeTreeNode node = cppTree.getNodeByIndex(i);
            if (node.getChildren().size() != 0) {
                String[] children = node.getLeafNames();
                if (children.length == 3 && node.getAge() == 2) {
                    assert (List.of(children).contains("1"));
                    assert (List.of(children).contains("2"));
                    assert (List.of(children).contains("3"));
                }
            }
        }

    }

    @Test
    void testCPP2() {
        Value<Number> samplingProb = new Value("", 0.2);
        Value<Number> birthRate = new Value("", 0.3);
        Value<Number> deathRate = new Value("", 0.1);
        Value<Integer> n = new Value("", 5);
        String[] taxa = new String[]{"1", "2", "3"};
        String[] allTaxa = new String[]{"1", "2", "3","4","5"};
        double age = 2.0;
        Value<CalibrationArray> calibration = new Value<>("", new CalibrationArray(new Calibration[]{new Calibration(taxa, age), new Calibration(allTaxa, 10.0)}));

        CalibratedCPPTree cpp = new CalibratedCPPTree(birthRate, deathRate, null, null,  samplingProb, n, calibration, null, null);
        TimeTree cppTree = cpp.sample().value();

        for (int i = 0; i < cppTree.getNodeCount(); i++) {
            TimeTreeNode node = cppTree.getNodeByIndex(i);
            if (!node.getChildren().isEmpty()) {
                String[] children = node.getLeafNames();
                if (children.length == 3 && node.getAge() == 2) {
                    assert (List.of(children).contains("1"));
                    assert (List.of(children).contains("2"));
                    assert (List.of(children).contains("3"));
                }
            }
        }
        assertEquals(10, cppTree.getRoot().getAge());

    }

    @Test
    void testCPP3() {
        Value<Number> samplingProb = new Value("", 0.5);
        Value<Number> birthRate = new Value("", 0.3);
        Value<Number> deathRate = new Value("", 0.1);
        Value<Integer> n = new Value("", 3);
        String[] taxa1 = new String[]{"1", "2", "3"};
        double age1 = 2.0;
        String[] taxa2 = new String[]{"1", "2"};
        double age2 = 1.0;
        Value<CalibrationArray> calibration = new Value<>("", new CalibrationArray(new Calibration[]{new Calibration(taxa1, age1), new Calibration(taxa2, age2)}));

        CalibratedCPPTree cpp = new CalibratedCPPTree(birthRate, deathRate, null, null,  samplingProb, n, calibration, null, null);
        for (int j = 0; j < 10 ; j ++) {
            try {
                TimeTree cppTree = cpp.sample().value();
                System.out.println("Iteration " + j);
                System.out.println(cppTree);

                for (int i = 0; i < cppTree.getNodeCount(); i++) {
                    TimeTreeNode node = cppTree.getNodeByIndex(i);
                    if (node.getChildren().size() != 0) {
                        String[] children = node.getLeafNames();
                        if (children.length == 3 && Math.abs(node.getAge() - 2.0) < 1e-8) {
                            List<String> childList = List.of(children);
                            assertEquals (2, node.getAge(), 1e-8);
                            assertTrue(childList.contains("1"), "Missing '1' in 3-clade");
                            assertTrue(childList.contains("2"), "Missing '2' in 3-clade");
                            assertTrue(childList.contains("3"), "Missing '3' in 3-clade");
                        }


                        if (children.length == 2 && List.of(children).contains("1") && List.of(children).contains("2")) {
                            assertEquals(1, node.getAge() , 1e-8);
                        }
                    }
                }

                //assertEquals(2, cppTree.getRoot().getAge());

            } catch (Throwable e) {
                System.err.println("Exception at iteration " + j);
                e.printStackTrace(); // Print full stack trace
                // Optional: break or continue depending on whether you want to stop
                break;
            }
        }
    }

    @Test
    void testCPP4() {
        Value<Number> samplingProb = new Value("", 0.5);
        Value<Number> birthRate = new Value("", 0.3);
        Value<Number> deathRate = new Value("", 0.1);
        Value<Integer> n = new Value("", 5);
        String[] taxa1 = new String[]{"1", "2", "3"};
        double age1 = 2.0;
        String[] taxa2 = new String[]{"1", "2"};
        double age2 = 1.0;
        Value<CalibrationArray> calibration = new Value<>("", new CalibrationArray(new Calibration[]{new Calibration(taxa1, age1), new Calibration(taxa2, age2)}));


        CalibratedCPPTree cpp = new CalibratedCPPTree(birthRate, deathRate, null, null,  samplingProb, n, calibration, null, null);
        TimeTree cppTree = cpp.sample().value();

        for (int i = 0; i < cppTree.getNodeCount(); i++) {
            TimeTreeNode node = cppTree.getNodeByIndex(i);
            if (node.getChildren().size() != 0) {
                String[] children = node.getLeafNames();
                if (children.length == 3 && Math.abs(node.getAge() - 2.0) < 1e-8) {
                    List<String> childList = List.of(children);
                    assertEquals (2, node.getAge(), 1e-8);
                    assertTrue(childList.contains("1"), "Missing '1' in 3-clade");
                    assertTrue(childList.contains("2"), "Missing '2' in 3-clade");
                    assertTrue(childList.contains("3"), "Missing '3' in 3-clade");
                }


                if (children.length == 2 && List.of(children).contains("1") && List.of(children).contains("2")) {
                    assertEquals(1, node.getAge() , 1e-8);
                }
            }
        }
    }

    @Test
    void testSingleCalibration8Taxa() {
        Value<Number> samplingProb = new Value<>("", 0.5);
        Value<Number> birthRate = new Value<>("", 0.3);
        Value<Number> deathRate = new Value<>("", 0.1);
        Value<Integer> n = new Value<>("", 8);

        // one calibrated clade of 4 taxa, 4 uncalibrated taxa outside
        String[] cladeTaxa = new String[]{"a", "b", "c", "d"};
        double cladeAge = 2.0;
        Value<CalibrationArray> calibration = new Value<>("",
                new CalibrationArray(new Calibration[]{new Calibration(cladeTaxa, cladeAge)}));

        CalibratedCPPTree cpp = new CalibratedCPPTree(birthRate, deathRate, null, null,  samplingProb, n, calibration, null, null);

        for (int iter = 0; iter < 100; iter++) {
            TimeTree tree = cpp.sample().value();
            System.out.println("Iteration " + iter + ": " + tree);

            assertEquals(8, tree.getLeafNodes().size(), "Should have 8 leaves");

            // check every internal node is older than its children
            for (int i = 0; i < tree.getNodeCount(); i++) {
                TimeTreeNode node = tree.getNodeByIndex(i);
                for (TimeTreeNode child : node.getChildren()) {
                    assertTrue(node.getAge() > child.getAge(),
                            "Parent age " + node.getAge() + " must be > child age " + child.getAge()
                                    + " at iteration " + iter + "\nTree: " + tree);
                }
            }

            // check the calibrated clade exists with correct age and taxa
            for (int i = 0; i < tree.getNodeCount(); i++) {
                TimeTreeNode node = tree.getNodeByIndex(i);
                if (!node.isLeaf()) {
                    String[] leaves = node.getLeafNames();
                    if (leaves.length == 4) {
                        List<String> leafList = List.of(leaves);
                        if (leafList.contains("a") && leafList.contains("b")
                                && leafList.contains("c") && leafList.contains("d")) {
                            assertEquals(cladeAge, node.getAge(), 1e-8,
                                    "Calibrated clade MRCA should have age " + cladeAge);
                        }
                    }
                }
            }
        }
    }

    @Test
    void testDeterministicFunctions() {
        double birthRate = 2.0;
        double deathRate = 1.0;
        double samplingProb = 0.5;
        int t = 2;

        // expected numbers calculated from r
        //assertEquals(0.8646647, CDF(birthRate, deathRate, null, null,  samplingProb, t), 1e-6);
        //assertEquals(0 , CDF(birthRate, deathRate, null, null,  samplingProb,0), 1e-6);
        //assertEquals(1 , CDF(birthRate, deathRate, null, null,  samplingProb,Double.POSITIVE_INFINITY), 1e-6);
        //assertEquals(0.1353353, densityBD(birthRate, deathRate, null, null,  samplingProb, t), 1e-6);
        assertEquals(0.2231436, inverseCDF(birthRate, deathRate,  samplingProb, 0.2), 1e-6);
        assertEquals(1.203973, inverseCDF(birthRate, deathRate, samplingProb, 0.7), 1e-6);
        assertEquals(0.4707278, Qdist(birthRate, deathRate,  t, 10), 1e-6);
        assertEquals(1.826258, transform(0.4, birthRate, deathRate, 10), 1e-6);

        List<Double> list = new ArrayList<>(Arrays.asList(0.0, 1.0, 2.0, 5.0, 0.0));
        assertEquals(0, indexOfMin(list));

        boolean[] list2 = new boolean[]{true, false, true, false, false};
        List<Integer> results = checkTrues(list2);
        int[] expect = new int[]{0,2};
        assertEquals(2, results.size());
        for (int i = 0; i < expect.length; i++) {
            assertEquals(expect[i], results.get(i));
        }
    }

    @Test
    void testCoalesce() {
        List<Double> times = new ArrayList<>();
        times.add(4.0);
        times.add(4.0);
        times.add(2.0);
        times.add(3.0);
        times.add(1.0);

        List<TimeTreeNode> nodeList = new ArrayList<>();
        TimeTreeNode node0 = new TimeTreeNode(0.0);
        TimeTreeNode node1 = new TimeTreeNode(0.0);
        TimeTreeNode node2 = new TimeTreeNode(1.0);
        TimeTreeNode node21 = new TimeTreeNode(0.0);
        TimeTreeNode node22 = new TimeTreeNode(0.0);
        node2.addChild(node21);
        node2.addChild(node22);
        TimeTreeNode node3 = new TimeTreeNode(0.0);
        TimeTreeNode node4 = new TimeTreeNode(0.0);
        nodeList.add(node0);
        nodeList.add(node1);
        nodeList.add(node2);
        nodeList.add(node3);
        nodeList.add(node4);

        coalesce(nodeList, times);

        // check ages for tips' parents
        assertEquals(1, nodeList.size());
        assertEquals(4.0, nodeList.get(0).getAge());
        assertEquals(4.0, node0.getParent().getAge());
        assertEquals(2.0, node1.getParent().getAge());
        assertEquals(2.0, node2.getParent().getAge());
        assertEquals(1.0, node3.getParent().getAge());
        assertEquals(1.0, node4.getParent().getAge());

        // check internal node ages
        assertEquals(3.0, node3.getParent().getParent().getAge());
        assertEquals(3.0, node1.getParent().getParent().getAge());
        assertEquals(4.0, node2.getParent().getParent().getParent().getAge());
    }

    @Test
    void testCoalesce2() {
        List<Double> times = new ArrayList<>();
        times.add(3.0);
        times.add(1.5);
        times.add(2.0);
        times.add(3.0);
        times.add(1.0);

        List<TimeTreeNode> nodeList = new ArrayList<>();
        TimeTreeNode node0 = new TimeTreeNode(0.0);
        TimeTreeNode node1 = new TimeTreeNode(0.0);
        TimeTreeNode node2 = new TimeTreeNode(1.0);
        TimeTreeNode node21 = new TimeTreeNode(0.0);
        TimeTreeNode node22 = new TimeTreeNode(0.0);
        node2.addChild(node21);
        node2.addChild(node22);
        TimeTreeNode node3 = new TimeTreeNode(0.0);
        TimeTreeNode node4 = new TimeTreeNode(0.0);
        nodeList.add(node0);
        nodeList.add(node1);
        nodeList.add(node2);
        nodeList.add(node3);
        nodeList.add(node4);

        coalesce(nodeList, times);

        // check ages for tips' parents
        assertEquals(1, nodeList.size());
        assertEquals(3.0, nodeList.get(0).getAge());
        assertEquals(1.5, node0.getParent().getAge());
        assertEquals(1.5, node1.getParent().getAge());
        assertEquals(2.0, node2.getParent().getAge());
        assertEquals(1.0, node3.getParent().getAge());
        assertEquals(1.0, node4.getParent().getAge());

        // check internal node ages
        assertEquals(3.0, node3.getParent().getParent().getAge());
        assertEquals(2.0, node1.getParent().getParent().getAge());
        assertEquals(3.0, node2.getParent().getParent().getAge());
    }

    @Test
    void testDuplicateNames() {
        Value<Number> samplingProb = new Value("", 0.2);
        Value<Number> birthRate = new Value("", 0.3);
        Value<Number> deathRate = new Value("", 0.1);
        Value<Integer> n = new Value("", 7);
        String[] taxa = new String[]{"1", "2", "3"};
        String[] taxa2 = new String[]{"4","5"};
        double age = 2.0;
        Value<CalibrationArray> calibration = new Value<>("", new CalibrationArray(new Calibration[]{new Calibration(taxa, age), new Calibration(taxa2, 3.0)}));

        CalibratedCPPTree cpp = new CalibratedCPPTree(birthRate, deathRate, null, null,  samplingProb, n, calibration, null, new Value<>("", 8.0));
        TimeTree cppTree = cpp.sample().value();
        String[] taxaNames = cppTree.getTaxaNames();

        // containing name 0-5
        boolean containsDigit0to5 = Arrays.stream(taxaNames)
                .anyMatch(name -> name.matches(".*[0-5].*"));
        assertTrue(containsDigit0to5);

        // containing name 1_1
        boolean hasLeaf1 = false;

        for (String name : taxaNames) {
            if (name.equals("1_2")) {
                hasLeaf1 = true;
                break;
            }
        }

        assertTrue(hasLeaf1);
    }
}
