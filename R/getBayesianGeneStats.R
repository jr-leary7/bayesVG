#' Fetch the table of per-gene HVG or SVG statistics.
#'
#' @name getBayesianGeneStats
#' @author Jack R. Leary
#' @description Given a user-provided object upon which the relevant HVG or SVG identification pipeline has been run, this function fetches the table of Bayesian gene statistics that is stored in the object's metadata.
#' @param obj An object of class \code{Seurat}, \code{SingleCellExperiment}, or \code{SpatialExperiment}. Defaults to NULL.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom SummarizedExperiment rowData
#' @importFrom Seurat DefaultAssay
#' @importFrom purrr reduce
#' @importFrom methods slot
#' @return A \code{data.frame} containing the relevant Bayesian gene statistics.
#' @export

getBayesianGeneStats <- function(obj = NULL) {
  # check inputs
  if (is.null(obj)) { cli::cli_abort("Please provide all inputs to getBayesianGeneStats().") }
  if (!(inherits(obj, "SingleCellExperiment") || inherits(obj, "Seurat") || inherits(obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class Seurat, SingleCellExperiment, or SpatialExperiment.") }
  # extract HVG or SVG summary table
  if (inherits(obj, "SingleCellExperiment")) {
    multi_subject_flag <- ifelse(is.null(obj@metadata$gene_stats_bayes), FALSE, TRUE)
    if (multi_subject_flag) {
      gene_summary <- purrr::reduce(obj@metadata$gene_stats_bayes, rbind)
    } else {
      gene_summary <- as.data.frame(SummarizedExperiment::rowData(obj))
    }
  } else if (inherits(obj, "SpatialExperiment")) {
    gene_summary <- as.data.frame(SummarizedExperiment::rowData(obj))
  } else if (inherits(obj, "Seurat")) {
    multi_subject_flag <- ifelse(is.null(obj@assays[[Seurat::DefaultAssay(obj)]]@misc$gene_stats_bayes), FALSE, TRUE)
    if (multi_subject_flag) {
      gene_summary <- purrr::reduce(obj@assays[[Seurat::DefaultAssay(obj)]]@misc$gene_stats_bayes, rbind)
    } else {
      version_check <- try({
        methods::slot(obj@assays[[Seurat::DefaultAssay(obj)]], name = "meta.data")
      }, silent = TRUE)
      if (inherits(version_check, "try-error")) {
        gene_summary <- obj@assays[[Seurat::DefaultAssay(obj)]]@meta.features
      } else {
        gene_summary <- obj@assays[[Seurat::DefaultAssay(obj)]]@meta.data
      }
    }
  }
  return(gene_summary)
}
