# R package imo: Ion Mirror Optimization
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

The R package provides functions to optimize the geometry and voltages of some planar or axially symmetric ion mirrors, with and without grids, where the axial potential can be calculated analytically.

## Installation
```r
if (!require("remotes")) { install.packages("remotes") }
remotes::install_github("pasturm/imo", auth_token = "github_pat")
```
where github_pat is a personal access token that can be generated in [https://github.com/settings/tokens?type=beta](https://github.com/settings/tokens?type=beta)

## Release notes
See the [NEWS file](https://github.com/pasturm/imo/blob/master/NEWS.md) for the latest release notes.

## Notes
* The tuning parameters are configured in a [configuration file](https://github.com/pasturm/imo/blob/master/inst/GLPM_config.toml).
* Design of experiments (DoE) and response surface methodology (RSM) is used for efficient optimization.
* The same approach is used in [SIMIONtuneR](https://github.com/pasturm/SIMIONtuneR) to efficiently optimize SIMION simulations.
