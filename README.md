
- [`bayesVG`](#bayesvg)
- [Installation](#installation)
- [Usage](#usage)
  - [Libraries](#libraries)
  - [HVG detection](#hvg-detection)
  - [SVG detection](#svg-detection)
- [Contact information](#contact-information)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `bayesVG`

<!-- badges: start -->

[![R-CMD-check](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml)
![last
commit](https://img.shields.io/github/last-commit/jr-leary7/bayesVG/main?color=darkgreen)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Coverage](https://codecov.io/gh/jr-leary7/bayesVG/graph/badge.svg)](https://app.codecov.io/gh/jr-leary7/bayesVG)
<!-- badges: end -->

# Installation

You can install the most recent version of `bayesVG` using:

``` r
remotes::install_github("jr-leary7/bayesVG")
```

# Usage

## Libraries

``` r
library(Seurat)
library(bayesVG)
```

## HVG detection

### Data

First, we load the 10X Genomics pbmc3k dataset, which is composed of
2,700 peripheral blood mononuclear cells from a single healthy donor.

``` r
data("seu_pbmc")
```

### Modeling

Now we’re able to model gene expression, summarize the posterior
distribution of variance for each gene, and classify the top 3000
most-variable genes as HVGs.

``` r
seu_pbmc <- findVariableFeaturesBayes(seu_pbmc, 
                                      n.cells.subsample = 1000L, 
                                      algorithm = "meanfield",
                                      save.model = TRUE) %>% 
            classifyHVGs(n.HVG = 3000L)
```

## SVG detection

### Data

First, we load the 10X Genomics anterior mouse brain dataset.

``` r
data("seu_brain")
```

Before running `bayesVG` for SVG detection it’s necessary to normalize
the expression data and identify a set of naive HVGs.

``` r
seu_brain <- SCTransform(seu_brain,
                         assay = "Spatial",
                         variable.features.n = 3000L,
                         vst.flavor = "v2",
                         return.only.var.genes = FALSE,
                         seed.use = 312,
                         verbose = FALSE)
```

### Modeling

Now we can model gene expression with an approximate Gaussian process,
summarize the spatial component of variance for each gene, and classify
the top 1000 most spatially variable genes as SVGs.

``` r
seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain, 
                                                kernel = "matern", 
                                                kernel.smoothness = 1.5, 
                                                algorithm = "meanfield", 
                                                n.cores = 4L, 
                                                save.model = TRUE) %>% 
             classifySVGs(n.SVG = 1000L)
```

# Contact information

This package is developed & maintained by Jack R. Leary. Feel free to
reach out by [opening an
issue](https://github.com/jr-leary7/bayesVG/issues) or by email
(<j.leary@ufl.edu>) if more detailed assistance is needed.
