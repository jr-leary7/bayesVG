% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichSpatialModules.R
\name{enrichSpatialModules}
\alias{enrichSpatialModules}
\title{Run enrichment analysis on SVG modules.}
\usage{
enrichSpatialModules(svg.clusters = NULL, species = "hsapiens")
}
\arguments{
\item{svg.clusters}{The results from \code{\link{clusterSVGsBayes}}. Defaults to NULL.}

\item{species}{A string specifying the species from which the cells originate. Defaults to "hsapiens".}
}
\value{
A \code{data.frame} with enrichment analysis results for each spatial module.
}
\description{
This function utilizes \code{\link[gprofiler2]{gost}} to perform enrichment analysis on each SVG module.
}
\examples{
data(seu_brain)
seu_brain <- Seurat::NormalizeData(seu_brain, verbose = FALSE) \%>\% 
             Seurat::FindVariableFeatures(nfeatures = 3000L, verbose = FALSE)
seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain,
                                                naive.hvgs = Seurat::VariableFeatures(seu_brain),
                                                kernel = "matern",
                                                kernel.smoothness = 1.5,
                                                algorithm = "meanfield",
                                                n.cores = 1L,
                                                save.model = TRUE) \%>\%
             classifySVGs(n.SVG = 1000L)
svg_clusters <- clusterSVGsBayes(seu_brain,
                                 svgs = Seurat::VariableFeatures(seu_brain),
                                 n.clust = 3L,
                                 n.cores = 1L)
enrich_res <- enrichSpatialModules(svg_clusters, species = "mmusculus")
}
\author{
Jack R. Leary
}
