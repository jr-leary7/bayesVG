#' Fetch the table of per-gene HVG or SVG statistics.
#'
#' @name getBayesianGeneStats
#' @author Jack R. Leary
#' @description Given a user-provided object upon which the relevant HVG or SVG identification pipeline has been run, this function fetches the table of Bayesian gene statistics that is stored in the object's metadata.
#' @param obj An object of class \code{Seurat}, \code{SingleCellExperiment}, or \code{SpatialExperiment}. Defaults to NULL.
#' @param sort.values A Boolean specifying whether the resulting \code{data.frame} should be sorted such that the most variable genes appear first. Defaults to TRUE.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom SummarizedExperiment rowData
#' @importFrom Seurat DefaultAssay
#' @importFrom purrr reduce
#' @importFrom methods slot
#' @importFrom dplyr filter arrange desc
#' @return A \code{data.frame} containing the relevant Bayesian gene statistics.
#' @export
#' @examples
#' data(seu_pbmc)
#' seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
#'                                       n.cells.subsample = 500L,
#'                                       algorithm = "meanfield",
#'                                       n.cores.per.chain = 1L,
#'                                       save.model = TRUE) %>% 
#'              classifyHVGs(n.HVG = 100L) 
#' gene_stats <- getBayesianGeneStats(seu_pbmc)

getBayesianGeneStats <- function(obj = NULL, sort.values = TRUE) {
  # check inputs
  if (is.null(obj)) { cli::cli_abort("Please provide all inputs to getBayesianGeneStats().") }
  if (!(inherits(obj, "SingleCellExperiment") || inherits(obj, "Seurat") || inherits(obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat}, {.pkg SingleCellExperiment}, or {.pkg SpatialExperiment}.") }
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
    multi_subject_flag <- ifelse(is.null(obj@assays[[Seurat::DefaultAssay(obj)]]@misc$gene_stats_bayes), 
                                 FALSE, 
                                 TRUE)
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
  spatial_flag <- ifelse("amplitude_mean" %in% colnames(gene_summary), 
                         TRUE, 
                         FALSE)
  if (spatial_flag) {
    gene_summary <- dplyr::filter(gene_summary, !is.na(amplitude_mean))
  }
  if (sort.values) {
    if (spatial_flag) {
      gene_summary <- dplyr::arrange(gene_summary, dplyr::desc(amplitude_mean))
    } else if ("dispersion_mean" %in% colnames(gene_summary)) {
      gene_summary <- dplyr::arrange(gene_summary, dplyr::desc(dispersion_mean))
    }
  }
  return(gene_summary)
}
