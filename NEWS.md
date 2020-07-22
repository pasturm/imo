# Version 0.1.0.9005

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
  

# Version 0.1.0

* Initial version.
