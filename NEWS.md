# bayesVG v0.0.7

+ Fixed a tiny bug in the Matern kernel code. 
+ Sped up QR decomposition of basis functions via conditional usage of LAPACK. 
+ Updated docs and made some functions less brittle. 
+ Updated some tests.
+ Updated the README. 
+ Added naive HVG fetching function. 
+ Made most of the Examples faster by decreasing the number of HVGs / SVGs. 
+ Removed a NOTE from R CMD check by declaring global variables used throughout the package.
+ Fixed some other minor bugs. 
+ Updated more tests. 
+ Added non-centered priors to all NB models. 
+ Removed dependency on `ggspavis`. 
+ Added support for plotting multiple genes at once to `plotSpatialExpression()`. 

# bayesVG v0.0.6 

+ Updated dependencies for R 4.5.1. 
+ Fixed failing tests related to `Seurat` to `SpatialExperiment` conversion. 
+ Updated spatial expression & attribute plotting functions to use `ggspavis::plotCoords()` as `ggspavis::plotSpots()` has been deprecated. 
+ Updated `plotModuleScores()` to account for changes to `ggplot2::geom_violin()` in `ggplot2` v4.0 with respect to how quantiles are drawn. 
+ Changed some dependencies in docs & made the code in Examples more performant.
+ Updated CITATION file to use new `bibentry()` format. 
+ Added non-centered prior parameterizations to the two Gaussian SVG models' Stan code. This increases accuracy as well as speed.
+ Dramatically sped up the SVG clustering model. 
+ Sped up matrix distance computation via the `fields` package and matrix scaling via the `coop` package. 
+ Made some of the Examples faster in order to speed up R CMD check. 

# bayesVG v0.0.5

+ Added function `plotSpatialExpression()` to generate a clean plot of gene expression overlaid on spatial coordinates given a `Seurat` or `SpatialExperiment` object. 
+ Added function `plotTissueImage()` to render a plot of the tissue image alone - with no information overlaid - given a `Seurat` object. 
+ Changed the way basis functions are generated so that they are orthonormal (mutually orthogonal and having unit norm) via a QR decomposition.
+ Added function `enrichSpatialModules()` to perform GSEA on clusters of SVGs using `gProfiler2` under the hood.
+ Made the default option to adjust for differing means among genes to be a hierarchical prior on the intercept instead of a fixed effect for gene depth. 
+ Added option to use the Negative-binomial likelihood for gene expression to `findSpatiallyVariableFeaturesBayes()`, along with accompanying Stan code. The relevant parameter is denoted `likelihood`. 
+ Updated test suite. 
+ Added function `plotModuleScores()` to generate spatial scatterplots, embedding scatterplots, and violin plots of spatial module scores.
+ Added alternative method of estimating global length-scale within `findSpatiallyVariableFeaturesBayes()` using aggregated variograms. 

# bayesVG v0.0.4 

+ Fixed bugs throughout. 
+ Added an R function and accompanying Stan code to soft-cluster SVG sets in a Bayesian manner. 
+ Upgraded clustering function and relevant Stan code to run in parallel and provide estimates of the log-likelihood and BIC of the final clustering. 
+ Updated documentation. 
+ Added multi-subject support to main HVG function and related downstream functions. 
+ Added new helper function `getBayesianGeneStats()` to pull HVG / SVG summary table from user-specified object's metadata.
+ Updated and expanded test suite. 
+ Updated `DESCRIPTION` file to point to BioConductor as this is needed for some dependencies e.g., `SpatialExperiment`.
+ Added function to compute spatial module scores for SVG clusters via the `UCell` package. 

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
