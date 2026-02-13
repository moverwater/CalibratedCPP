# Likelihood Benchmark and Validation

We benchmark the computation time for the likelihood calculation of a single fully balanced tree with 128 leaves as we increase the number of cherries for which we want to condition the age. The time is compared with the [Heled and Drummond (2015)](https://doi.org/10.1093/sysbio/syu089) implementation of the 'Calibrated Birth-Death process'.

We also validate the likelihood computation by comparing log-likelihood values for a fixed tree varying the birth-rate with both the Heled and Drummond and Calibrated CPP likelihood implementations.

## Results

**Likelihood Benchmark: CPP v Heled and Drummond**
<!-- start table -->

|   Number of Calibrations |   Time (μs) Calibrated CPP |   Time (μs) Heled & Drummond |
|--------------------------|----------------------------|------------------------------|
|                      1.0 |                       2.47 |                         6.45 |
|                      2.0 |                       3.51 |                       227.66 |
|                      3.0 |                       6.07 |                      9580.50 |
|                      4.0 |                      14.95 |                    342198.91 |
|                      5.0 |                      47.06 |                  11440218.12 |
|                      6.0 |                     153.94 |                 267042755.72 |

<!-- end table -->

**Benchmark and Validation Plot**
![](./combined_benchmark_validation.png)

## Usage

Run benchmark (takes approximately 50 minutes):
```bash
mvn clean test-compile
java -Xss32m -cp "$(mvn dependency:build-classpath | grep -v '\[INFO\]' | tr '\n' ':'):target/classes:target/test-classes" \
    calibratedcpp.LikelihoodBenchmark
```


 Run validation:
```bash
mvn test-compile
mvn test -Dtest=CalibratedCoalescentPointProcessTest#heledAndDrummondComparison
```

 Generate plot:
 ```bash
python validation/calibratedcpp/validation_and_benchmark/plot_output.py
 ```

## Output Files

 - [`benchmark_results.csv`](./benchmark_results.csv).
 - [`validation_results.csv`](./validation_results.csv).
 - [`combined_benchmark_validation.png`](./combined_benchmark_validation.png).
