# bayesVG v0.0.6 

+ Updated dependencies for R 4.5.1. 
+ Fixed failing tests related to `Seurat` to `SpatialExperiment` conversion. 
+ Updated spatial expression & attribute plotting functions to use `ggspavis::plotCoords()` as `ggspavis::plotSpots()` has been deprecated. 
+ Updated `plotModuleScores()` to account for changes to `ggplot2::geom_violin()` with respect to how quantiles are drawn. 

# bayesVG v0.0.5

+ Added function `plotSpatialExpression()` to generate a clean plot of gene expression overlaid on spatial coordinates given a `Seurat` or `SpatialExperiment` object. 
+ Added function `plotTissueImage()` to render a plot of the tissue image alone - with no information overlaid - given a `Seurat` object. 
+ Changed the way basis functions are generated so that they are orthonormal (mutually orthogonal and having unit norm) via a QR decomposition.
+ Added function `enrichSpatialModules()` to perform GSEA on clusters of SVGs using `gProfiler2` under the hood.
+ Made the default option to adjust for differing means among genes to be a hierarchical prior on the intercept instead of a fixed effect for gene depth. 
+ Added option to use the Negative-binomial likelihood for gene expression to `findSpatiallyvariableFeaturesBayes()`, along with accompanying Stan code. The relevant parameter is denoted `likelihood`. 
+ Updated test suite. 
+ Added function `plotModuleScores()` to generate spatial scatterplots, embedding scatterplots, and violin plots of spatial module scores.
+ Added alternative method of estimating global length-scale within `findSpatiallyvariableFeaturesBayes()` using aggregated variograms. 

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
