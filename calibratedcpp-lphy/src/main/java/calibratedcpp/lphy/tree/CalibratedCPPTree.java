package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TaxaConditionedTreeGenerator;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.*;
import static calibratedcpp.lphy.tree.CPPUtils.*;

public class CalibratedCPPTree extends TaxaConditionedTreeGenerator implements GenerativeDistribution<TimeTree> {

    Value<Number> birthRate;
    Value<Number> deathRate;
    Value<CalibrationArray> calibrations;
    Value<Number> stemAge;
    Value<Number> rho;
    Value<String[]> otherNames;
    double conditionAge;
    boolean rootConditioned = false;

    List<String> nameList;
    public static final String calibrationsName = "calibrations";
    public static final String stemAgeName = "stemAge";
    public static final String otherTaxaNames = "otherNames";

    public CalibratedCPPTree(@ParameterInfo(name = BirthDeathConstants.lambdaParamName, description = "per-lineage birth rate.") Value<Number> birthRate,
                             @ParameterInfo(name = BirthDeathConstants.muParamName, description = "per-lineage death rate.") Value<Number> deathRate,
                             @ParameterInfo(name = BirthDeathConstants.rhoParamName, description = "sampling probability") Value<Number> rho,
                             @ParameterInfo(name = DistributionConstants.nParamName, description = "the total number of taxa.") Value<Integer> n,
                             @ParameterInfo(name = calibrationsName, description = "an array of calibrations generated from a MRCA prior (i.e. ConditionedMRCAPrior or MRCAPrior)") Value<CalibrationArray> calibrations,
                             @ParameterInfo(name = otherTaxaNames, description = "a string array of taxa names for non-calibrated tips", optional = true) Value<String[]> otherNames,
                             @ParameterInfo(name = stemAgeName, description = "the stem age working as condition time", optional = true) Value<Number> stemAge) {
        super(n, null, null);
        // check legal params
        if (calibrations == null){
            throw new NullPointerException("Calibrations should not be null!");
        }
        if (stemAge != null) {
            LoggerUtils.log.warning("Stem age will be ignored if root calibration is provided.");
        }

        this.rho = rho;
        this.calibrations = calibrations;
        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.otherNames = otherNames;
        this.stemAge = stemAge;
    }

    @GeneratorInfo(name = "CalibratedCPP", examples = {"CalibratedCPPTree.lphy"},
            description = "The Calibrated Coalescent Point Process (calibrated CPP) method accepts one or more clade taxa and generates a tip-labelled time tree. If a root age is provided, the method is conditioned on root age. If the stem age is provided, the origin is the stem age.")
    @Override
    public RandomVariable<TimeTree> sample() {
        // obtain pass in parameters
        double birthRate = getBirthRate().value().doubleValue();
        double deathRate = getDeathRate().value().doubleValue();
        double samplingProb = getSamplingProb().value().doubleValue();
        int n = getN().value();
        CalibrationArray calibrationArray = getCalibrations().value();
        Calibration[] calibrations = calibrationArray.getCalibrationArray();

        // initialise params
        double rootAge = 0;
        TimeTree tree = new TimeTree();

        List<String> backUpNames = new ArrayList<>();

        // step1: get valid clade calibrations
        List<Calibration> cladeCalibrations = new ArrayList<>(Arrays.stream(calibrations).toList());

        // sort it with decreasing order
        cladeCalibrations.sort((c1, c2) -> Double.compare(c2.getAge(), c1.getAge()));

        // if root calibration is already in clade calibration
        if (cladeCalibrations.get(0).getTaxa().length == n){
            rootConditioned = true;
            rootAge = cladeCalibrations.get(0).getAge();
            // if only one root calibration, then return cpp
            if (cladeCalibrations.size() == 1){
                CPPTree cpp = new CPPTree(getBirthRate(), getDeathRate(), getSamplingProb(),
                        new Value<>("", cladeCalibrations.get(0).getTaxa()), getN(), new Value<>("", cladeCalibrations.get(0).getAge()), null);
                tree = cpp.sample().value();
                return new RandomVariable<>("", tree, this);
            } else {
                // else remove the root calibration from cladeCalibrations
                backUpNames.addAll(Arrays.asList(cladeCalibrations.get(0).getTaxa()));
                cladeCalibrations.remove(cladeCalibrations.get(0));
            }
        }

        // step2: get all maximal calibration
        List<Calibration> maximalCalibrations = getMaximalCalibrations(cladeCalibrations);

        // map the taxa names for calibration clades
        int index = 0;
        int cladeSizes = 0;
        String[][] taxaNames = new String[maximalCalibrations.size()][];
        for (Calibration entry : maximalCalibrations) {
            taxaNames[index] = entry.getTaxa();
            cladeSizes += taxaNames[index].length;
            index++;
        }

        // calculate the number of nodes
        // m = non-clade tips + clade roots
        int m = n - cladeSizes + maximalCalibrations.size();

        /* initialise the lists
            A : holding the indices of inactive nodes (wait for assign)
            l : holding the indices of active nodes (has assigned)
            times : holding the times of internal nodes, the first element is root or stem age
            nodeAges : holding the times for each node
            nodeList : holding all nodes
         */
        List<Integer> A = new ArrayList<>(m);
        for (int i = 0; i < m; i++) {
            A.add(i);
        }
        int[] l = new int[m];
        List<Double> times = new ArrayList<>(Collections.nCopies(m, 0.0));
        List<Double> nodeAges = new ArrayList<>(Collections.nCopies(m, 0.0));
        List<TimeTreeNode> nodeList = new ArrayList<>((Collections.nCopies(m, null)));

        // step3: calculate condition age (root or stem age)
        // if rootConditioned, then condition on root
        // if !rootConditioned, then use stem age or sample one
        if (rootConditioned) {
            int ind;
            ind = random.nextInt(m - 1) + 1; // [1, m-1]
            if (m == 2){
                ind = 1;
            }

            times.set(ind, rootAge);
            conditionAge = rootAge;
        } else {
            if (getStemAge()!= null) {
                conditionAge = getStemAge().value().doubleValue();
            } else {
                conditionAge = simRandomStem(birthRate, deathRate, maximalCalibrations.get(0).getAge(), n);
            }
        }

        // step4: build clades for each maximalCalibrations
        nameList = new ArrayList<>(n);

        // set name list
        // TODO: check if name provided overlap the automatic names make them have prefix
        for (int i = 0; i< maximalCalibrations.size(); i++) {
            String[] currentNames = maximalCalibrations.get(i).getTaxa();
            String[] uniqueNames = new String[currentNames.length];
            for (int j = 0; j< currentNames.length; j++) {
                //String newName = "clade" + i + "_" + maximalCalibrationsEntries.get(i).getValue()[j];
                String newName = currentNames[j];

                nameList.add(newName);
                uniqueNames[j] = newName;
                backUpNames.remove(uniqueNames[j]);
            }
            maximalCalibrations.get(i).setTaxa(uniqueNames);
        }

        // loop through all maximalCalibrations
        for (int i = 0; i < maximalCalibrations.size(); i++) {
            // step1: get subclades
            List<Calibration> subClades = getNestedClades(maximalCalibrations.get(i), cladeCalibrations);
            // step2: get sampled element
            // calculate weights
            double w = CDF(birthRate, deathRate, samplingProb, conditionAge) -
                    CDF(birthRate, deathRate, samplingProb, cladeCalibrations.get(i).getAge());
            // calculate score s for each node
            int[] s = calculateScore(A, m, times);
            // calculate weight for each node
            double[] weights = getWeights(s, w);

            if (A.size() == 1){
                l[i] = A.get(0);
            } else {
                // sample one element from A with probability weights
                l[i] = A.get(sampleIndex(weights));
            }

            // step3: construct subtrees
            // construct calibration clade leaf names
            Calibration[] clades = new Calibration[subClades.size()];

            for (int j = 0; j < subClades.size(); j++) {
                Calibration cal = new Calibration(subClades.get(j).getTaxa());
                cal.setAge(subClades.get(j).getAge());
                clades[j] = cal;
            }

            // simulate a tree for these clades, only offer calibrations
            CalibratedCPPTree calibratedCPPTree = new CalibratedCPPTree(getBirthRate(),
                    getDeathRate(), getSamplingProb(),
                    new Value<>("n", maximalCalibrations.get(i).getTaxa().length),
                    new Value<>("", new CalibrationArray(clades)), null, null);

            // put clade mrca into a list waiting for assign
            TimeTree subTree = calibratedCPPTree.sample().value();
            nodeList.set(l[i], subTree.getRoot());

            // step4: assign unresolved node times for l[i]
            // once done, remove l[i] from list A
            // deal with the nodes have l[i] still 0
            if (times.get(l[i]) == 0) {
                double time = sampleTimes(birthRate, deathRate, samplingProb, maximalCalibrations.get(i).getAge(), conditionAge, 1)[0];
                times.set(l[i],time);
            }

            if (l[i] < m -1  && times.get(l[i] + 1) == 0 ) {
                double time = sampleTimes(birthRate, deathRate, samplingProb, maximalCalibrations.get(i).getAge(), conditionAge, 1)[0];
                times.set(l[i] + 1, time);
            }

            nodeAges.set(l[i], maximalCalibrations.get(i).getAge());

            // remove corresponding node in A
            A.remove(Integer.valueOf(l[i]));
        }

        // step5: organise times
        // after calibrations, sample times for remaining unassigned nodes
        List<Integer> zeroIndices = new ArrayList<>();
        for (int i = 0; i < times.size(); i++) {
            if (times.get(i) == 0) zeroIndices.add(i);
        }
        double[] sampledTimes = sampleTimes(birthRate, deathRate, samplingProb, 0, conditionAge, zeroIndices.size());
        for (int j = 0; j < zeroIndices.size(); j++) {
            times.set(zeroIndices.get(j), sampledTimes[j]);
        }

        // set the first node to be the max, make it the root
        if (times.size() > 1) {
            double max = times.get(1);
            for (int i = 1; i < times.size() ; i++) {
                if (times.get(i) > max) {
                    max = times.get(i);
                }
            }
            times.set(0, max); // set this the largest
        }

        // or if not root conditioned, set it to stem age
        if (!rootConditioned){
            times.set(0, conditionAge);
        }

        // step6: fill in nodelist
        // get non-clade taxa
        List<String> outGroupTaxa = new ArrayList<>();
        if (getOtherNames() != null){
            String[] otherNames = getOtherNames().value();
            for (String otherName : otherNames) {
                nameList.add(otherName);
                outGroupTaxa.add(otherName);
            }
        } else {
            if (!backUpNames.isEmpty()){
                for (String name : backUpNames) {
                    nameList.add(name);
                    outGroupTaxa.add(name);
                }
            }

        }

        // fit other names in if there are non-assigned names
        int nameListSize = nameList.size();
        for (int i = 0; i < n - nameListSize; i++) {
            nameList.add(String.valueOf(i));
            outGroupTaxa.add(String.valueOf(i));
        }

        // get random order for non-clade taxa
        Collections.shuffle(outGroupTaxa);

        // Assign remaining uncalibrated taxa names to available node positions
        int ind = 0;
        for (int i = 0; i < nodeList.size() && ind < outGroupTaxa.size(); i++) {
            if (nodeList.get(i) == null) {
                TimeTreeNode tip = new TimeTreeNode(0);
                tip.setId(outGroupTaxa.get(ind));
                nodeList.set(i, tip);
                ind++;
            }
        }

        // step7: coalesce
        // combine sub-CPPs into the final tree
        while (nodeList.size() > 1) {
            // start from the youngest node
            int j = indexOfMin(times);
            if (times.size() == 2) {
                j = 1; // make it being the second one if there are only 2 nodes
            }

            // build relationship
            TimeTreeNode child_left = nodeList.get(j-1);
            TimeTreeNode child_right = nodeList.get(j);

            TimeTreeNode parent = new TimeTreeNode(times.get(j));
            parent.addChild(child_left);
            parent.addChild(child_right);

            child_left.setParent(parent);
            child_right.setParent(parent);

            // give ages
            child_left.setAge(nodeAges.get(j-1));
            child_right.setAge(nodeAges.get(j));

            // adjust indices
            nodeList.set(j-1, parent);
            nodeList.remove(j);

            // set parent node to time of current node
            nodeAges.set(j-1, times.get(j));

            // remove the time and age of the second node
            times.remove(j);
            nodeAges.remove(j);
        }

        tree.setRoot(nodeList.get(0), true);
        return new RandomVariable<>("CPPTree", tree, this);
    }

    /*
       Functions
    */
    private static int[] calculateScore(List<Integer> A, int m, List<Double> times) {
        int[] s = new int[A.size()];
        for (int j = 0; j < A.size(); j++) {
            int nodeIndex = A.get(j);
            int count = 0;
            // Check if i < m and i+1 is within bounds
            if (nodeIndex < m - 1 && times.get(nodeIndex + 1) == 0) {
                count++;
            }
            // Check if nodeIndex is within bounds and times[i] == 0
            if (times.get(nodeIndex) == 0) {
                count++;
            }
            s[j] = count;
        }
        return s;
    }


    private static double[] getWeights(int[] s, double w) {
        double sumOfWeights = 0;
        double[] weights = new double[s.length];

        for (int i = 0; i < s.length; i++) {
            weights[i] = Math.pow(w, s[i]);
            sumOfWeights += weights[i];
        }

        // normalise weights
        for (int i = 0; i < s.length; i++) {
            weights[i] /= sumOfWeights;
        }

        return weights;
    }


    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        map.put(BirthDeathConstants.lambdaParamName, birthRate);
        map.put(BirthDeathConstants.muParamName, deathRate);
        map.put(BirthDeathConstants.rhoParamName, rho);
        map.put(DistributionConstants.nParamName, n);
        map.put(calibrationsName, calibrations);
        if (stemAge != null) map.put(stemAgeName, stemAge);
        if (otherNames != null) map.put(otherTaxaNames, otherNames);
        return map;
    }

    public void setParam(String paramName, Value value){
        if (paramName.equals(BirthDeathConstants.lambdaParamName)) birthRate = value;
        else if (paramName.equals(BirthDeathConstants.muParamName)) deathRate = value;
        else if (paramName.equals(BirthDeathConstants.rhoParamName)) rho = value;
        else if (paramName.equals(DistributionConstants.nParamName)) n = value;
        else if (paramName.equals(calibrationsName)) calibrations = value;
        else if (paramName.equals(otherTaxaNames)) otherNames = value;
        else if (paramName.equals(stemAgeName)) stemAge = value;
        else throw new IllegalArgumentException("Unknown parameter name: " + paramName);
    }

    public Value<Integer> getN(){
        return getParams().get(DistributionConstants.nParamName);
    }

    public Value<Number> getBirthRate(){
        return getParams().get(BirthDeathConstants.lambdaParamName);
    }

    public Value<Number> getDeathRate(){
        return getParams().get(BirthDeathConstants.muParamName);
    }

    public Value<Number> getSamplingProb(){
        return getParams().get(BirthDeathConstants.rhoParamName);
    }

    public Value<CalibrationArray> getCalibrations(){
        return getParams().get(calibrationsName);
    }

    public Value<Number> getStemAge(){return getParams().get(stemAgeName);}

    public Value<String[]> getOtherNames(){
        return getParams().get(otherTaxaNames);
    }

    public Value<Double> getOrigin(){
        return new Value<>("", conditionAge);
    }

    public boolean getRootCondition(){
        return rootConditioned;
    }
}

