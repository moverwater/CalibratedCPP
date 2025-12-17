# Likelihood validation comparison

We compute likelihood values for a fixed tree with 100 leaves and set of 4 calibrations (including the root) using both the [Heled and Drummond (2015)](https://doi.org/10.1093/sysbio/syu089) 'Calibrated Birth-Death model' and the 'Calibrated CPP'. To run compare the likelihood values:
 - Go to the [src/testjava/calibratedcpp/](../../../src/test/java/calibratedcpp) directory in the calibratedcpp-beast subproject.
 - Run the `heledAndDrummondComparison()` method in `CalibratedCoalescentPointProcessTest.java`.
 - This will write a csv file to the current directory called [comparison_results.csv](comparison_results.csv)
