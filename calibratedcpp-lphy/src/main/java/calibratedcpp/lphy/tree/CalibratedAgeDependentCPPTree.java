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

public class CalibratedAgeDependentCPPTree extends TaxaConditionedTreeGenerator implements GenerativeDistribution<TimeTree> {

    Value<Number> birthRate;
    Value<Number> rho;
    Value<Number> lifetime;
    Value<CalibrationArray> calibrations;
    Value<Number> stemAge;
    Value<String[]> otherNames;
    double conditionAge;
    boolean rootConditioned = false;

    List<String> nameList;
    Set<String> usedNames;

    public static final String calibrationsName = "calibrations";
    public static final String stemAgeName     = "stemAge";
    public static final String otherTaxaNames  = "otherNames";
    public static final String lifetimeName    = "lifetime";

    public CalibratedAgeDependentCPPTree(
            @ParameterInfo(name = BirthDeathConstants.lambdaParamName,
                    description = "per-lineage birth rate.") Value<Number> birthRate,
            @ParameterInfo(name = BirthDeathConstants.rhoParamName,
                    description = "sampling probability.") Value<Number> rho,
            @ParameterInfo(name = DistributionConstants.nParamName,
                    description = "the total number of taxa.") Value<Integer> n,
            @ParameterInfo(name = lifetimeName,
                    description = "individual lifetime distribution; must come from a continuous distribution " +
                            "(e.g. Gamma, LogNormal, Exp). The sampled value is used as a point estimate of " +
                            "the mean lifetime for prior-predictive tree simulation. The full distribution is " +
                            "passed to BEAST for exact likelihood evaluation.") Value<Number> lifetime,
            @ParameterInfo(name = calibrationsName,
                    description = "an array of calibrations generated from a MRCA prior " +
                            "(i.e. ConditionedMRCAPrior or MRCAPrior).") Value<CalibrationArray> calibrations,
            @ParameterInfo(name = otherTaxaNames,
                    description = "a string array of taxa names for non-calibrated tips.",
                    optional = true) Value<String[]> otherNames,
            @ParameterInfo(name = stemAgeName,
                    description = "the stem age working as condition time.",
                    optional = true) Value<Number> stemAge) {
        super(n, null, null);

        if (calibrations == null) {
            throw new NullPointerException("Calibrations should not be null!");
        }
        if (lifetime == null) {
            throw new NullPointerException("lifetime should not be null!");
        }

        if (stemAge != null && calibrations.value().getCalibrationArray()[0].getTaxa().length == n.value()) {
            LoggerUtils.log.warning("Stem age will be ignored if root calibration is provided.");
        }

        this.birthRate   = birthRate;
        this.rho         = rho;
        this.lifetime    = lifetime;
        this.calibrations = calibrations;
        this.otherNames  = otherNames;
        this.stemAge     = stemAge;
    }

    @GeneratorInfo(name = "CalibratedAgeDependentCPP", examples = {},
            description = "The Calibrated Coalescent Point Process with age-dependent individual lifetimes. "
                    + "The lifetime parameter must come from a continuous distribution (Gamma, LogNormal, or Exp). "
                    + "Tree sampling uses effectiveDeathRate = 1/lifetime for initialisation; "
                    + "the exact age-dependent likelihood is evaluated in BEAST via CalibratedAgeDependentBirthDeathModel.")
    @Override
    public RandomVariable<TimeTree> sample() {
        double birthRateVal = getBirthRate().value().doubleValue();
        // 1/mean_lifetime as an effective death rate for BD-based tree initialisation
        double effectiveDeathRate = 1.0 / getLifetime().value().doubleValue();
        double samplingProb = getSamplingProb().value().doubleValue();
        int n = getN().value();
        CalibrationArray calibrationArray = getCalibrations().value();
        Calibration[] calibrations = calibrationArray.getCalibrationArray();

        double rootAge = 0;
        TimeTree tree = new TimeTree();
        List<String> backUpNames = new ArrayList<>();

        // step1: sorted clade calibrations, descending by age
        List<Calibration> cladeCalibrations = new ArrayList<>(Arrays.stream(calibrations).toList());
        cladeCalibrations.sort((c1, c2) -> Double.compare(c2.getAge(), c1.getAge()));

        if (cladeCalibrations.getFirst().getTaxa().length == n) {
            rootConditioned = true;
            rootAge = cladeCalibrations.getFirst().getAge();
            if (cladeCalibrations.size() == 1) {
                CPPTree cpp = new CPPTree(
                        new Value<>("", birthRateVal), new Value<>("", effectiveDeathRate),
                        null, null, getSamplingProb(),
                        new Value<>("", cladeCalibrations.getFirst().getTaxa()), getN(),
                        new Value<>("", cladeCalibrations.getFirst().getAge()), null);
                tree = cpp.sample().value();
                return new RandomVariable<>("", tree, this);
            } else {
                backUpNames.addAll(Arrays.asList(cladeCalibrations.getFirst().getTaxa()));
                cladeCalibrations.remove(cladeCalibrations.getFirst());
            }
        }

        // step2: maximal calibrations
        List<Calibration> maximalCalibrations = getMaximalCalibrations(cladeCalibrations);

        int index = 0, cladeSizes = 0;
        String[][] taxaNames = new String[maximalCalibrations.size()][];
        for (Calibration entry : maximalCalibrations) {
            taxaNames[index] = entry.getTaxa();
            cladeSizes += taxaNames[index].length;
            index++;
        }
        int m = n - cladeSizes + maximalCalibrations.size();

        List<Integer> A = new ArrayList<>(m);
        for (int i = 0; i < m; i++) A.add(i);
        int[] l = new int[m];
        List<Double> times    = new ArrayList<>(Collections.nCopies(m, 0.0));
        List<TimeTreeNode> nodeList = new ArrayList<>(Collections.nCopies(m, null));

        // step3: condition age
        if (rootConditioned) {
            int ind = (m == 2) ? 1 : random.nextInt(m - 1) + 1;
            times.set(ind, rootAge);
            conditionAge = rootAge;
        } else {
            if (getStemAge() != null) {
                conditionAge = getStemAge().value().doubleValue();
            } else {
                int idx = 0;
                while (Double.isNaN(conditionAge) || Double.isInfinite(conditionAge) || conditionAge == 0.0) {
                    conditionAge = simRandomStem(birthRateVal, effectiveDeathRate,
                            maximalCalibrations.getFirst().getAge(), n);
                    if (++idx > 200) {
                        throw new RuntimeException(
                                "The stem age cannot be sampled with the current parameters. Please provide a stemAge.");
                    }
                }
            }
        }

        // step4: build clades
        nameList = new ArrayList<>(n);
        usedNames = new HashSet<>();

        for (Calibration maximalCalibration : maximalCalibrations) {
            String[] currentNames = maximalCalibration.getTaxa();
            String[] uniqueNames  = new String[currentNames.length];
            for (int j = 0; j < currentNames.length; j++) {
                String newName = makeUnique(currentNames[j], usedNames);
                nameList.add(newName);
                uniqueNames[j] = newName;
                backUpNames.remove(currentNames[j]);
            }
            maximalCalibration.setTaxa(uniqueNames);
        }

        for (int i = 0; i < maximalCalibrations.size(); i++) {
            List<Calibration> subClades = getNestedClades(maximalCalibrations.get(i), cladeCalibrations);

            double w = CDF(birthRateVal, effectiveDeathRate, samplingProb, conditionAge)
                    - CDF(birthRateVal, effectiveDeathRate, samplingProb, maximalCalibrations.get(i).getAge());
            int[]    s       = calculateScore(A, m, times);
            double[] weights = getWeights(s, w);

            l[i] = (A.size() == 1) ? A.getFirst() : A.get(sampleIndex(weights));

            double   cladeAge    = maximalCalibrations.get(i).getAge();
            String[] subcladeTaxa = maximalCalibrations.get(i).getTaxa();

            Calibration[] clades = new Calibration[subClades.size() + 1];
            clades[0] = maximalCalibrations.get(i);
            for (int j = 0; j < subClades.size(); j++) {
                Calibration cal = new Calibration(subClades.get(j).getTaxa());
                cal.setAge(subClades.get(j).getAge());
                clades[j + 1] = cal;
            }

            // recursive: subclade also uses the same lifetime distribution
            CalibratedAgeDependentCPPTree subTreeGen = new CalibratedAgeDependentCPPTree(
                    new Value<>("", birthRateVal), getSamplingProb(),
                    new Value<>("n", subcladeTaxa.length),
                    getLifetime(),
                    new Value<>("", new CalibrationArray(clades)),
                    null, null);
            nodeList.set(l[i], subTreeGen.sample().value().getRoot());

            if (l[i] > 0 && times.get(l[i]) == 0) {
                times.set(l[i],
                        sampleTimes(birthRateVal, effectiveDeathRate, samplingProb, cladeAge, conditionAge, 1)[0]);
            }
            if (l[i] < m - 1 && times.get(l[i] + 1) == 0) {
                times.set(l[i] + 1,
                        sampleTimes(birthRateVal, effectiveDeathRate, samplingProb, cladeAge, conditionAge, 1)[0]);
            }
            A.remove(Integer.valueOf(l[i]));
        }

        // step5: fill remaining times
        List<Integer> zeroIndices = new ArrayList<>();
        for (int i = 1; i < times.size(); i++) {
            if (times.get(i) == 0) zeroIndices.add(i);
        }
        double[] sampledTimes = sampleTimes(birthRateVal, effectiveDeathRate, samplingProb,
                0, conditionAge, zeroIndices.size());
        for (int j = 0; j < zeroIndices.size(); j++) {
            times.set(zeroIndices.get(j), sampledTimes[j]);
        }

        if (times.size() > 2) {
            double max = times.get(1);
            for (int i = 1; i < times.size(); i++) if (times.get(i) > max) max = times.get(i);
            times.set(0, max);
        } else if (times.size() == 2) {
            times.set(0, times.get(1));
        } else {
            throw new RuntimeException("Unreachable");
        }

        if (rootConditioned && Math.abs(times.getFirst() - rootAge) > 1e-8) {
            throw new RuntimeException("Max age does not match root age when root-conditioned.");
        }
        if (!rootConditioned) {
            times.set(0, conditionAge);
        }

        // step6: non-clade taxa
        List<String> outGroupTaxa = new ArrayList<>();
        if (getOtherNames() != null) {
            for (String name : getOtherNames().value()) {
                String newName = makeUnique(name, usedNames);
                nameList.add(newName);
                outGroupTaxa.add(newName);
            }
        } else {
            for (String name : backUpNames) {
                nameList.add(name);
                outGroupTaxa.add(name);
            }
        }
        int nameListSize = nameList.size();
        for (int i = 0; i < n - nameListSize; i++) {
            String newName = makeUnique(String.valueOf(i), usedNames);
            nameList.add(newName);
            outGroupTaxa.add(newName);
        }
        Collections.shuffle(outGroupTaxa);

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
        CalibratedCPPTree.coalesce(nodeList, times);
        tree.setRoot(nodeList.getFirst(), true);
        if (!rootConditioned && conditionAge != nodeList.getFirst().getAge()) {
            tree.getRoot().setRootStem(times.getFirst());
        }

        return new RandomVariable<>("AgeDependentCPPTree", tree, this);
    }

    private static int[] calculateScore(List<Integer> A, int m, List<Double> times) {
        int[] s = new int[A.size()];
        for (int j = 0; j < A.size(); j++) {
            int idx = A.get(j);
            int count = 0;
            if (idx < m - 1 && times.get(idx + 1) == 0) count++;
            if (idx > 0     && times.get(idx)     == 0) count++;
            s[j] = count;
        }
        return s;
    }

    private static double[] getWeights(int[] s, double w) {
        double sum = 0;
        double[] weights = new double[s.length];
        for (int i = 0; i < s.length; i++) { weights[i] = Math.pow(w, s[i]); sum += weights[i]; }
        for (int i = 0; i < s.length; i++) weights[i] /= sum;
        return weights;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        map.put(BirthDeathConstants.lambdaParamName, birthRate);
        map.put(BirthDeathConstants.rhoParamName, rho);
        map.put(lifetimeName, lifetime);
        map.put(calibrationsName, calibrations);
        if (stemAge   != null) map.put(stemAgeName,    stemAge);
        if (otherNames != null) map.put(otherTaxaNames, otherNames);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if      (paramName.equals(BirthDeathConstants.lambdaParamName))   birthRate    = value;
        else if (paramName.equals(BirthDeathConstants.rhoParamName))       rho          = value;
        else if (paramName.equals(lifetimeName))                            lifetime     = value;
        else if (paramName.equals(DistributionConstants.nParamName))        n            = value;
        else if (paramName.equals(calibrationsName))                        calibrations = value;
        else if (paramName.equals(otherTaxaNames))                          otherNames   = value;
        else if (paramName.equals(stemAgeName))                             stemAge      = value;
        else throw new IllegalArgumentException("Unknown parameter name: " + paramName);
    }

    public Value<Integer>         getN()           { return getParams().get(DistributionConstants.nParamName); }
    public Value<Number>          getBirthRate()   { return getParams().get(BirthDeathConstants.lambdaParamName); }
    public Value<Number>          getSamplingProb(){ return getParams().get(BirthDeathConstants.rhoParamName); }
    public Value<Number>          getLifetime()    { return getParams().get(lifetimeName); }
    public Value<CalibrationArray>getCalibrations(){ return getParams().get(calibrationsName); }
    public Value<Number>          getStemAge()     { return getParams().get(stemAgeName); }
    public Value<String[]>        getOtherNames()  { return getParams().get(otherTaxaNames); }
    public Value<Double>          getOrigin()      { return new Value<>("", conditionAge); }
    public boolean                getRootCondition(){ return rootConditioned; }
}
