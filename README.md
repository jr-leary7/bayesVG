
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
[![CodeFactor](https://www.codefactor.io/repository/github/jr-leary7/bayesvg/badge/main)](https://www.codefactor.io/repository/github/jr-leary7/bayesvg/overview/main)
<!-- badges: end -->

# Installation

You can install the most recent version of `bayesVG` using:

``` r
remotes::install_github("jr-leary7/bayesVG")
```

# Usage

## Libraries

``` r
library(dplyr)
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
                                      n.cells.subsample = 500L, 
                                      algorithm = "meanfield",
                                      save.model = TRUE) %>% 
            classifyHVGs(n.HVG = 3000L)
```

We can extract the summary table and classify the top 3,000 genes as
HVGs like so. These genes can then be used as the basis for downstream
analyses such as PCA, clustering, UMAP visualization, etc.

``` r
summary_hvg <- arrange(seu_pbmc@assays$RNA@meta.data, desc(dispersion_mean))
top3k_hvgs <- summary_hvg$gene[1:3000]
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

Now we can model gene expression with an approximate multivariate
hierarchical Gaussian process (GP), summarize the spatial component of
variance for each gene, and classify the top 1000 most spatially
variable genes as SVGs.

``` r
seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain, 
                                                naive.hvgs = VariableFeatures(seu_brain), 
                                                kernel = "matern", 
                                                kernel.smoothness = 1.5, 
                                                algorithm = "meanfield", 
                                                n.cores = 4L, 
                                                save.model = TRUE) %>% 
             classifySVGs(n.SVG = 1000L)
```

We can extract the summary table and classify the top 1,000 genes as
SVGs like so. These genes can then be used as the basis for downstream
analyses such as PCA, clustering, UMAP visualization, etc.

``` r
summary_svg <- arrange(seu_brain@assays$SCT@meta.features, amplitude_mean_rank)
top1k_svgs <- summary_svg$gene[1:1000]
```

# Contact information

This package is developed & maintained by Jack R. Leary. Feel free to
reach out by [opening an
issue](https://github.com/jr-leary7/bayesVG/issues) or by email
(<j.leary@ufl.edu>) if more detailed assistance is needed.
