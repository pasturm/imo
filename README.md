# R package imo: Ion Mirror Optimization
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

The R package provides functions to optimize the geometry and voltages of planar or axially symmetric ion mirrors, where the potential can be calculated analytically. Design of expriments (DoE) and response surface methodology (RSM) is used for efficient optimization. 


## Installation
``` r
if (!require("devtools")) { install.packages("devtools") }
devtools::install_github("pasturm/imo")
```

## Release notes
See the [NEWS file](https://github.com/pasturm/imo/blob/master/NEWS.md) for the latest release notes.


## Notes
* Parallel computing is used to speed up the optimization.
* The tuning parameters are configured in a [configuration file](https://github.com/pasturm/imo/blob/master/inst/GLPMtuneR_config.toml)
* The design of experiments and response surface method follows the approach 
of the TOFWERK Thuner and underlying MKS MODDE-Q software. Notable differences to Thuner/MODDE are:
    * It is open source (+). 
    * The response surface model optimization works much better (due to improved optimization algorithms and desirability functions) (+).
    * The optimizaition is much easier to configure (+). 
    * It can be used to efficiently optimize ion mirrors where the potential can be calculated analytically (+).
    * It can be used to efficiently optimize SIMION simulations ([SIMIONtuneR](https://github.com/pasturm/SIMIONtuneR)) (+).
    * It is not a self-contained program and does not have a graphical user interface (-). 
