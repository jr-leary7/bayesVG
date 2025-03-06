#' Cluster SVGs into modules.
#'
#' @name clusterSVGsBayes
#' @author Jack R. Leary
#' @description This downstream analysis function clusters identified SVGs into modules via an approximate Bayesian soft-clustering approach.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} upon which \code{\link{findSpatiallyVariableFeaturesBayes}} and \code{\link{classifySVGs}} have been run. Defaults to NULL.
#' @param svgs A character vector containing the identified SVGs. Defaults to NULL. 
#' @param n.clust An integer specifying the number of clusters to fit to the data. Defaults to 5. 
#' @import magrittr
#' @importFrom SingleCellExperiment logcounts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @return A \code{data.frame} containing the per-SVG soft cluster assignments.
#' @export

clusterSVGsBayes <- function(sp.obj = NULL, 
                             svgs = NULL, 
                             n.clust = 5L) {
  # check inputs 
  if (is.null(sp.obj) || is.null(svgs)) { stop("All arguments to clusterSVGsBayes() must be supplied.") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { stop("Please provide an object of class Seurat or SpatialExperiment.") }
  # extract matrix of normalized gene expression
  if (inherits(sp.obj, "Seurat")) {
    expr_mtx <- Seurat::GetAssayData(sp.obj,
                                     assay = Seurat::DefaultAssay(sp.obj),
                                     layer = "data")
  } else {
    expr_mtx <- SingleCellExperiment::logcounts(sp.obj)
  }
  expr_mtx <- expr_mtx[svgs, ]
}
