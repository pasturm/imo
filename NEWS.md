# Version 0.1.1

* `run_imo(tuneR_config, type = "GLPM")` normalizes the focal distance to the 
  height H.
* Added `write` parameter in `run_imo()` which (if `FALSE`) allows to 
  run the optimization without writing output files to disk.
* Improved plotting of results.
* Removed the `plot_coeffs()` and `plot_results()` functions.
* `run_imo()` returns the optimized best point of the last run.
* Updated the configuration file templates.
* All distances in GLPM functions are normalized to the height H of the mirror
  electrodes.
* The GLPM optimization reports the resolution obtained with 10 % energy 
  variation.
* Added `plot` parameter in `run_imo()` which (if `FALSE`) does not plot 
  anything.
* Takes starting values from `bestpoint` variable if `resume=TRUE` and no 
  "bestpoint_run.txt" file is found.
* Added `digits` parameter in `run_imo()` which controls the number
  of decimal places to print when printing the best point values.
* Removed parallelization as it is unstable.
* Fixed bug where not all electrodes were considered.

# Version 0.1.0

* Initial version.
