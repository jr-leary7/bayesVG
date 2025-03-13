# bayesVG v0.0.4 

+ Fixed bugs throughout. 
+ Added an R function and accompanying Stan code to soft-cluster SVG sets in a Bayesian manner. 
+ Upgraded clustering function and relevant Stan code to run in parallel and provide estimates of the log-likelihood and BIC of the final clustering. 
+ Updated documentation. 
+ Added multi-subject support to main HVG function and related downstream functions. 
+ Added new helper function `getBayesianGeneStats()` to pull HVG / SVG summary table from user-specified object's metadata.
+ Updated and expanded test suite. 
+ Updated `DESCRIPTION` file to point to BioConductor as this is needed for some dependencies e.g., `SpatialExperiment`.

# bayesVG v0.0.3 

+ Added function to plot HVG results. 
+ Fixed bugs throughout. 
+ Expanded test suite. 
+ Added option (and corresponding Stan code) to SVG modeling function to adjust for total gene sequencing depth. 
+ Switched to usage of the `cli` R package for verbose messages instead of the default `message()` function. 

# bayesVG v0.0.2 

+ Finished primary development of the HVG and SVG modeling functions. 
+ Updated and improved (with respect to efficiency and accuracy) the Stan code used in the SVG modeling function. 
+ Added a README. 
+ Improved function documentation. 
+ Changed structure of Stan code and compiled executable. 
+ Squashed some bugs. 
+ Added some (very) basic tests.
+ Started adding GitHub Actions for CI/CD. 
+ Added support for `SpatialExperiment` objects to SVG modeling function.

# bayesVG v0.0.1

+ Initial package structure. 
+ Implemented functions to identify HVGs in scRNA-seq data and SVGs in spatial transcriptomics data. 
+ Added `NEWS.md` file to track changes. 
