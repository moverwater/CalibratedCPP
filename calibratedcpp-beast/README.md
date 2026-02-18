# Calibrated Coalescent Point Process BEAST2 Implementation

Birth-death models are typically parameterised by a per-lineage birth-rate $(\lambda)$, a per-lineage death-rate $(\mu)$, and a probability with which each individual at the present time is sampled $(\rho)$. Other parameters of interest are:
- The reproductive number: $R=\lambda/\mu$.
- The diversification rate: $d=\lambda-\mu$.
- Turnover: $\tau=\mu/\lambda$.

We can specify any two of these parameters together except for turnover $(\tau)$ and the reproductive number $(R)$.

Both `BirthDeathModel` and `BirthDeathSkylineModel` extend `CalibratedCoalescentPointProcess` which itself extends `SpeciesTreeDistribution`. The  inputs of a `CalibratedCoalescentPointProcess` object are:
- `origin`: A `RealParameter` for the start time of the process. This is optional.
- `conditionOnRoot`: A `Boolean` which is true the process is conditioned on the root height. Either `origin` must be provided or `conditionOnRoot` should be set to true.
- `tree`: A phylogenetic tree.
- `calibrations`: A list of `CalibrationClade` objects which consists of a `TaxonSet` on which we want to condition being monophyletic and the TMRCA of the clade. Note: providing the all the taxa in the tree as a `calibration` is equivalent to setting `conditionOnRoot=true`.
- `conditionOnCalibrations`: A `Boolean` which is `true` if we want to condition on the `calibrations` provided and `false` otherwise. Default is `true`.

## `BirthDeathModel`

The birth-death model takes as inputs any TWO of birth-rate $(\lambda)$, death-rate $(\mu)$, diversification rate $(d)$, reproductive number $(R)$, and turnover $(\tau)$ as `RealParameter` objects; and the sampling probability $(\rho)$ as a `RealParameter` in the interval $[0,1]$.

## `BirthDeathSkylineModel`

The birth-death skyline model takes as inputs any TWO of the time dependent piecewise constant parameters: birth-rate $(\lambda(t))$, death-rate $(\mu(t))$, diversification rate $(d(t))$, effective reproductive number $(R_e(t))$, and turnover $(\tau(t))$ as `SkylineParameter` objects; and the sampling probability $(\rho)$ as a `RealParameter` in the interval $[0,1]$.

`SkylineParameter` objects take as inputs:
- Rates as `RealParameter` objects where rates are ordered from root to tips.
- Optionally, change times can be specified as `RealParameter`objects. If this is not provided, the rates are spread evenly over the age of the tree (origin/root height).
- Reverse is a `Boolean` which is `true` if the change times are given as ages before the present. The default is `false`.  
- Relative is a `Boolean` which is `true` if the change times are given as relative to the age of the tree (origin/root height). The default is `false`.