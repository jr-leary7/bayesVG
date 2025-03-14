% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getBayesianGeneStats.R
\name{getBayesianGeneStats}
\alias{getBayesianGeneStats}
\title{Fetch the table of per-gene HVG or SVG statistics.}
\usage{
getBayesianGeneStats(obj = NULL, sort.values = TRUE)
}
\arguments{
\item{obj}{An object of class \code{Seurat}, \code{SingleCellExperiment}, or \code{SpatialExperiment}. Defaults to NULL.}

\item{sort.values}{A Boolean specifying whether the resulting \code{data.frame} should be sorted such that the most variable genes appear first. Defaults to TRUE.}
}
\value{
A \code{data.frame} containing the relevant Bayesian gene statistics.
}
\description{
Given a user-provided object upon which the relevant HVG or SVG identification pipeline has been run, this function fetches the table of Bayesian gene statistics that is stored in the object's metadata.
}
\author{
Jack R. Leary
}
