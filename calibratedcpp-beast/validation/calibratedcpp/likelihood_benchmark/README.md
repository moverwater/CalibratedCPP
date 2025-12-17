# Likelihood computation benchmark

We benchmark the computation time for the likelihood calculation of a single fully balanced tree with 128 leaves as we increase the number of cherries for which we want to condition the age. The time is compared with the [Heled and Drummond (2015)](https://doi.org/10.1093/sysbio/syu089) implementation of the `Calibrated Birth-Death process'.

To run the benchmark comparison:
 - Go to the [src/test/java/calibratedcpp/ directory](../../../src/test/java/calibratedcpp/) in the calibratedcpp-beast sub-project.
 - Run `LikelihoodBenchmark.java`.
 - This will write a csv file called [benchmark_results.csv](./benchmark_results.csv) to this directory.
