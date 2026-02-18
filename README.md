# The Calibrated Coalescent Point Process
[![Docker CI](https://github.com/moverwater/CalibratedCPP/actions/workflows/maven-tests.yml/badge.svg)](https://github.com/moverwater/CalibratedCPP/actions/workflows/maven-tests.yml)

This is an implementation of the Calibrated Coalescent Point Process (Calibrated CPP) in BEAST2.

The CPP is a model of ultrametric trees where node ages are i.i.d. random variables. This captures a general class of birth-death processes with time-dependent birth and death rates (as in BDSKY) and age dependent death rates.
- Individuals die at some rate $\mu(x,t)$ that depends on the age of the individual $x$ and the time $t$.
- Individuals give birth to new individuals at a rate $\lambda(t)$ depending only on time.

[Lambert and Stadler (2012)](https://doi.org/10.1016/j.tpb.2013.10.002) show that these models give a uniform distribution over ranked labelled (or oriented) tree topologies, and that the node ages are i.i.d. random variables. The node age density where $\lambda=\lambda(t)$ and $\mu=\mu(x,t)$ and each individual is sampled with probability $\rho$ is,

$$q(t) = \frac{\rho\lambda (\lambda-\mu)e^{-(\lambda-\mu)t}}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}$$

the cumulative distribution function is,

$$Q(t) = \frac{\rho\lambda(1-e^{-(\lambda-\mu)t})}{\rho\lambda+(\lambda(1-\rho)-\mu)e^{-(\lambda-\mu)t}}.$$

The Calibrated CPP is a calibrated tree prior using the CPP. Calibrated tree priors are used for molecular clock dating by conditioning on the existence and ages of the most recent common ancestors of monophyletic clades.

## Project structure
CalibratedCPP contains 3 subproject:

### calibratedcpp-beast

The implementation has the following structure:
- The abstract class `CoalescentPointProcessModel` has abstract methods `calculateNodeAgeDensity()` and `calculateNodeAgeCDF()` for the density and CDF of the node age. This class extends `SpeciesTreeDistribution` and takes a list of `calibrations`, and the `origin` age OR `conditionOnRoot` as inputs.
- `BirthDeathModel` extends `CoalescentPointProcessModel` and implements `calculateNodeAgeDensity()` and `calculateNodeAgeCDF()` with node age density and CDF for the constant rate birth-death process.
- `BirthDeathSkylineModel` extends `CoalescentPointProcessModel` and implements `calculateNodeAgeDensity()` and `calculateNodeAgeCDF()` with node age density and CDF for the birth-death process with piecewise constant rates.
- `CalibratedCoalescentPointProcess` extends `SpeciesTreeDistribution` and takes a list of `calibrations`, and the `origin` age OR `conditionOnRoot` as inputs.

### calibratedcpp-lphy

The Lphy simulator is currently implemented within LinguaPhylo as a generative method called `CalibratedCPP`:
- The simulator takes birth rate, death rate, sampling probability for the birth-death process.
- Calibration information is taken from the output of `ConditionedMRCAPrior`, other leaf names are optional to pass in.
- Stem age can be passed in if there is no root calibration, but if root calibration and stemAge are both specified, the tree will still be root conditioned.

### calibratedcpp-lphybeast

This subproject is converting Lphy simulators to XMLs for BEAST2 running:
- The simulator set conditionOnCalibrations as true for default, users should manually modify this flag in the output XML to make it turn off.
- Construct BirthDeathModel with birth death parameters and conditions on the origin of the tree.
- Set calibrations in CalibrationPrior, taking upper and lower bounds of the calibration nodes that passed to ConditionedMRCAPrior.

## License

CalibratedCPP is free software.  It is distributed under the terms of version 3 of the GNU General Public License.  A copy of this license should be found in the file [COPYING](./COPYING) located in the root directory of this repository. If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.













