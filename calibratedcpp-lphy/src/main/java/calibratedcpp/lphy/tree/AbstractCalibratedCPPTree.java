package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TaxaConditionedTreeGenerator;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;

import java.util.*;

import static calibratedcpp.lphy.tree.CPPUtils.*;

/**
 * Shared calibrated-CPP tree-building logic (clade decomposition, weighted slot
 * assignment, coalescing) for a coalescent point process whose node ages come
 * from some birth-death law. Concrete subclasses only need to supply that law
 * via {@link #cdf}, {@link #sampleAges} and {@link #sampleStemAge} (mirroring
 * {@code calculateLogNodeAgeCDF}/{@code calculateLogNodeAgeDensity} on the BEAST
 * side, see {@code CalibratedCoalescentPointProcess}), plus a constant-rate pair
 * for the no-calibration fallback that delegates to plain {@link CPPTree}.
 */
public abstract class AbstractCalibratedCPPTree extends TaxaConditionedTreeGenerator implements GenerativeDistribution<TimeTree> {

    Value<Number> rho;
    Value<CalibrationArray> calibrations;
    Value<String[]> otherNames;
    Value<Number> stemAge;
    Value<Number> rootAge;
    double conditionAge;
    boolean rootConditioned = false;

    List<String> nameList;
    Set<String> usedNames;

    public static final String calibrationsName = "calibrations";
    public static final String stemAgeName = "stemAge";
    public static final String otherTaxaNames = "otherNames";
    public static final String rootAgeName = "rootAge";

    protected AbstractCalibratedCPPTree(Value<Integer> n, Value<Number> rho, Value<CalibrationArray> calibrations,
                                         Value<String[]> otherNames, Value<Number> stemAge, Value<Number> rootAge) {
        super(n, null, null);

        if (calibrations == null && rootAge == null) {
            throw new IllegalArgumentException("Either calibrations or rootAge must be provided.");
        }

        this.rho = rho;
        this.calibrations = calibrations;
        this.otherNames = otherNames;
        this.stemAge = stemAge;
        this.rootAge = rootAge;
    }

    // ****** node age density hooks — the only thing concrete subclasses implement ******

    /** Recompute any resolved rate values from the current parameter bindings. Called once per sample(). */
    protected abstract void resolveRates();

    /** CDF of the node age law up to time t, given the currently resolved rates. */
    protected abstract double cdf(double t);

    /** Draw nSims node ages from the law truncated to [lowerTime, upperTime]. */
    protected abstract double[] sampleAges(double lowerTime, double upperTime, int nSims);

    /** Sample a random stem age older than greaterThan for nTaxa taxa. */
    protected abstract double sampleStemAge(double greaterThan, int nTaxa);

    /** Birth/death rate pair to hand to the plain constant-rate {@link CPPTree} fallback. */
    protected abstract Value<Number> getConstantBirthRateValue();

    protected abstract Value<Number> getConstantDeathRateValue();

    /** Construct a subclade generator of the same concrete type, reusing this generator's rate parameters. */
    protected abstract AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations);

    private RandomVariable<TimeTree> delegateToCPPTree(Value<String[]> taxa, Value<Number> conditionRootAge) {
        CPPTree cpp = new CPPTree(getConstantBirthRateValue(), getConstantDeathRateValue(), null, null,
                getSamplingProb(), taxa, getN(), conditionRootAge, null);
        return new RandomVariable<>("", cpp.sample().value(), this);
    }

    @Override
    public RandomVariable<TimeTree> sample() {
        resolveRates();

        // if no calibrations, delegate to CPPTree conditioned on rootAge
        if (getCalibrations() == null) {
            return delegateToCPPTree(getOtherNames(), getRootAge());
        }

        // obtain pass in parameters
        double samplingProb = getSamplingProb().value().doubleValue();
        int n = getN().value();
        CalibrationArray calibrationArray = getCalibrations().value();
        Calibration[] calibrations = calibrationArray.getCalibrationArray();

        // initialise params
        double rootAgeValue = 0;
        TimeTree tree = new TimeTree();

        List<String> backUpNames = new ArrayList<>();

        // step1: get valid clade calibrations
        List<Calibration> cladeCalibrations = new ArrayList<>(Arrays.stream(calibrations).toList());

        // sort it with decreasing order
        cladeCalibrations.sort((c1, c2) -> Double.compare(c2.getAge(), c1.getAge()));

        // if root calibration is already in clade calibration
        if (cladeCalibrations.getFirst().getTaxa().length == n) {
            rootConditioned = true;
            rootAgeValue = cladeCalibrations.getFirst().getAge();
            // if only one root calibration, then return cpp
            if (cladeCalibrations.size() == 1) {
                return delegateToCPPTree(new Value<>("", cladeCalibrations.getFirst().getTaxa()),
                        new Value<>("", cladeCalibrations.getFirst().getAge()));
            } else {
                // else remove the root calibration from cladeCalibrations
                backUpNames.addAll(Arrays.asList(cladeCalibrations.getFirst().getTaxa()));
                cladeCalibrations.remove(cladeCalibrations.getFirst());
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
            nodeList : holding all nodes
         */
        List<Integer> A = new ArrayList<>(m);
        for (int i = 0; i < m; i++) {
            A.add(i);
        }
        int[] l = new int[m];
        List<Double> times = new ArrayList<>(Collections.nCopies(m, 0.0));
        List<TimeTreeNode> nodeList = new ArrayList<>((Collections.nCopies(m, null)));

        // step3: calculate condition age (root or stem age)
        // if rootConditioned, then condition on root
        // if !rootConditioned, then use stem age or sample one
        if (rootConditioned) {
            int ind;
            ind = random.nextInt(m - 1) + 1; // [1, m-1]
            if (m == 2) {
                ind = 1;
            }

            times.set(ind, rootAgeValue);
            conditionAge = rootAgeValue;
        } else {
            if (getStemAge() != null) {
                conditionAge = getStemAge().value().doubleValue();
            } else {
                int idx = 0;
                while (Double.isNaN(conditionAge) || Double.isInfinite(conditionAge) || conditionAge == 0.0) {
                    conditionAge = sampleStemAge(maximalCalibrations.getFirst().getAge(), n);
                    idx++;
                    if (idx > 200) {
                        throw new RuntimeException("The stem age cannot be sampled because of the bad parameter combination. Please provide a stemAge.");
                    }
                }
            }
        }

        // step4: build clades for each maximalCalibrations
        nameList = new ArrayList<>(n);
        usedNames = new HashSet<>();

        for (Calibration maximalCalibration : maximalCalibrations) {
            String[] currentNames = maximalCalibration.getTaxa();
            String[] uniqueNames = new String[currentNames.length];
            for (int j = 0; j < currentNames.length; j++) {
                String newName = makeUnique(currentNames[j], usedNames);

                nameList.add(newName);
                uniqueNames[j] = newName;
                backUpNames.remove(currentNames[j]);
            }
            maximalCalibration.setTaxa(uniqueNames);
        }

        // loop through all maximalCalibrations
        for (int i = 0; i < maximalCalibrations.size(); i++) {
            // step1: get subclades
            List<Calibration> subClades = getNestedClades(maximalCalibrations.get(i), cladeCalibrations);
            // step2: get sampled element
            // calculate weights
            double w = cdf(conditionAge) - cdf(maximalCalibrations.get(i).getAge());
            // calculate score s for each node
            int[] s = calculateScore(A, m, times);
            // calculate weight for each node
            double[] weights = getWeights(s, w);

            if (A.size() == 1) {
                l[i] = A.getFirst();
            } else {
                // sample one element from A with probability weights
                l[i] = A.get(sampleIndex(weights));
            }

            // step3: construct subtrees
            double cladeAge = maximalCalibrations.get(i).getAge();
            String[] subcladeTaxa = maximalCalibrations.get(i).getTaxa();

            Calibration[] clades = new Calibration[subClades.size() + 1];
            clades[0] = maximalCalibrations.get(i);
            for (int j = 0; j < subClades.size(); j++) {
                Calibration cal = new Calibration(subClades.get(j).getTaxa());
                cal.setAge(subClades.get(j).getAge());
                clades[j + 1] = cal;
            }

            // simulate a tree for these clades, only offer calibrations
            AbstractCalibratedCPPTree subTreeGen = newSubClade(subcladeTaxa.length, new CalibrationArray(clades));

            // put clade into nodeList
            TimeTree subTree = subTreeGen.sample().value();
            nodeList.set(l[i], subTree.getRoot());

            // sample times at l[i] and l[i]+1 conditioned to be older than the clade age,
            if (l[i] > 0 && times.get(l[i]) == 0) {
                times.set(l[i], sampleAges(cladeAge, conditionAge, 1)[0]);
            }

            if (l[i] < m - 1 && times.get(l[i] + 1) == 0) {
                times.set(l[i] + 1, sampleAges(cladeAge, conditionAge, 1)[0]);
            }

            // remove corresponding node in A
            A.remove(Integer.valueOf(l[i]));
        }

        // step5: organise times
        // after calibrations, sample times for remaining unassigned nodes
        List<Integer> zeroIndices = new ArrayList<>();
        for (int i = 1; i < times.size(); i++) {
            if (times.get(i) == 0) zeroIndices.add(i);
        }
        double[] sampledTimes = sampleAges(0, conditionAge, zeroIndices.size());
        for (int j = 0; j < zeroIndices.size(); j++) {
            times.set(zeroIndices.get(j), sampledTimes[j]);
        }

        // set the first node to be the max, make it the root
        if (times.size() > 2) {
            double max = times.get(1);
            for (int i = 1; i < times.size(); i++) {
                if (times.get(i) > max) {
                    max = times.get(i);
                }
            }
            times.set(0, max); // set this the largest
        } else if (times.size() == 2) {
            times.set(0, times.get(1));
        } else {
            throw new RuntimeException("Unreachable");
        }

        if (rootConditioned) {
            if (Math.abs(times.getFirst() - rootAgeValue) > 1e-8) {
                // shouldn't have this thrown theoretically
                throw new RuntimeException("The max age is not root age when root conditioned");
            }
        }

        // or if not root conditioned, set it to stem age
        if (!rootConditioned) {
            times.set(0, conditionAge);
        }

        // step6: fill in nodelist
        // get non-clade taxa
        List<String> outGroupTaxa = new ArrayList<>();
        if (getOtherNames() != null) {
            String[] otherNamesArr = getOtherNames().value();
            for (String otherName : otherNamesArr) {
                String newName = makeUnique(otherName, usedNames);
                nameList.add(newName);
                outGroupTaxa.add(newName);
            }
        } else {
            if (!backUpNames.isEmpty()) {
                for (String name : backUpNames) {
                    nameList.add(name);
                    outGroupTaxa.add(name);
                }
            }
        }

        // fit other names in if there are non-assigned names
        int nameListSize = nameList.size();
        for (int i = 0; i < n - nameListSize; i++) {
            String indexName = String.valueOf(i);
            String newName = makeUnique(indexName, usedNames);
            nameList.add(newName);
            outGroupTaxa.add(newName);
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
        coalesce(nodeList, times);

        tree.setRoot(nodeList.getFirst(), true);
        if (!rootConditioned && conditionAge != nodeList.getFirst().getAge()) {
            tree.getRoot().setRootStem(times.getFirst());
        }

        return new RandomVariable<>("CPPTree", tree, this);
    }

    // public for unit test
    public static void coalesce(List<TimeTreeNode> nodeList, List<Double> times) {
        while (nodeList.size() > 1) {
            // start from the youngest node
            int j = indexOfMin(times);
            if (times.size() == 2) {
                j = 1; // make it being the second one if there are only 2 nodes
            }

            // build relationship
            TimeTreeNode child_left = nodeList.get(j - 1);
            TimeTreeNode child_right = nodeList.get(j);

            TimeTreeNode parent = new TimeTreeNode(times.get(j));
            parent.addChild(child_left);
            parent.addChild(child_right);

            child_left.setParent(parent);
            child_right.setParent(parent);

            // adjust indices
            nodeList.set(j - 1, parent);
            nodeList.remove(j);

            // remove the time and age of the second node
            times.remove(j);
        }
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
            if (nodeIndex < m - 1 && times.get(nodeIndex + 1) == 0) count++;
            // Check if nodeIndex is within bounds and times[i] == 0
            if (nodeIndex > 0 && times.get(nodeIndex) == 0) count++;
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
        map.put(BirthDeathConstants.rhoParamName, rho);
        map.put(DistributionConstants.nParamName, n);
        if (calibrations != null) map.put(calibrationsName, calibrations);
        if (rootAge != null) map.put(rootAgeName, rootAge);
        if (stemAge != null) map.put(stemAgeName, stemAge);
        if (otherNames != null) map.put(otherTaxaNames, otherNames);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(BirthDeathConstants.rhoParamName)) rho = value;
        else if (paramName.equals(DistributionConstants.nParamName)) n = value;
        else if (paramName.equals(calibrationsName)) calibrations = value;
        else if (paramName.equals(otherTaxaNames)) otherNames = value;
        else if (paramName.equals(stemAgeName)) stemAge = value;
        else if (paramName.equals(rootAgeName)) rootAge = value;
        else throw new IllegalArgumentException("Unknown parameter name: " + paramName);
    }

    public Value<Integer> getN() {
        return getParams().get(DistributionConstants.nParamName);
    }

    public Value<Number> getSamplingProb() {
        return getParams().get(BirthDeathConstants.rhoParamName);
    }

    public Value<CalibrationArray> getCalibrations() {
        return getParams().get(calibrationsName);
    }

    public Value<Number> getStemAge() {
        return getParams().get(stemAgeName);
    }

    public Value<Number> getRootAge() {
        return getParams().get(rootAgeName);
    }

    public Value<String[]> getOtherNames() {
        return getParams().get(otherTaxaNames);
    }

    public Value<Double> getOrigin() {
        return new Value<>("", conditionAge);
    }

    public boolean getRootCondition() {
        return rootConditioned;
    }
}
