#' Classify genes as spatially variable.
#'
#' @name classifySVGs
#' @author Jack R. Leary
#' @description After estimating per-gene spatial variation statistics, classify genes as spatially variable or not in a variety of ways.
#' @param sp.obj An object of class \code{Seurat} containing spatial data. Defaults to NULL.
#' @param selection.method A string specifying what method should be used to classify genes as SVGs. Must be one of "rank", "quantile", or "cutoff". Defaults to "rank".
#' @param n.SVG An integer specifying the number of SVGs to select (if using rank-based selection). Defaults to 1000.
#' @param quantile.SVG A double specifying the quantile cutoff used to classify SVGs (if using quantile-based selection). Defaults to 0.75.
#' @param cutoff A double specifying the cutoff value for spatial variation (depending on how \code{selection.variable} is defined) used to classify SVGs (if using cutoff-based selection). Defaults to 0.05.
#' @import magrittr
#' @importFrom Seurat DefaultAssay VariableFeatures
#' @importFrom dplyr select arrange desc slice_head pull mutate filter if_else
#' @importFrom stats quantile
#' @return An object of class \code{Seurat} with SVG metadata added.
#' @seealso \code{\link{findSpatiallyVariableFeaturesBayes}}
#' @seealso \code{\link[SeuratObject]{SVFInfo}}
#' @export

classifySVGs <- function(sp.obj = NULL,
                         selection.method = "rank",
                         n.SVG = 1000L,
                         quantile.SVG = 0.75,
                         cutoff = 0.05) {
  # check inputs
  if (is.null(sp.obj)) { stop("Please provide an object to classifySVGs().") }
  selection.method <- tolower(selection.method)
  if (!selection.method %in% c("rank", "quantile", "cutoff")) { stop("Please provide a valid SVG selection method to classifyHVGs().") }
  # extract gene spatial variation statistics
  version_check <- try({
    slot(sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]], name = "meta.data")
  }, silent = TRUE)
  if (inherits(version_check, "try-error")) {
    gene_summary <- sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.features
  } else {
    gene_summary <- sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.data
  }
  # identify HVGs based on user-specified method
  if (selection.method == "rank") {
    svgs <- dplyr::arrange(gene_summary, mean_rank) %>%
            dplyr::slice_head(n = n.HVG) %>%
            dplyr::pull(gene)
  } else if (selection.method == "quantile") {
    quantile_cutoff <- as.numeric(stats::quantile(gene_summary[, "mean"], quantile.HVG))
    svgs <- dplyr::arrange(gene_summary, mean_rank) %>%
            dplyr::filter(mean >= quantile_cutoff) %>%
            dplyr::pull(gene)
  } else if (selection.method == "cutoff") {
    svgs <- dplyr::arrange(gene_summary, mean_rank) %>%
            dplyr::filter(mean >= cutoff) %>%
            dplyr::pull(gene)
  }
  # add HVG classification back to object metadata
  gene_summary <- dplyr::mutate(gene_summary, svg = dplyr::if_else(gene %in% svgs, TRUE, FALSE))
  if (inherits(version_check, "try-error")) {
    sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.features <- gene_summary
  } else {
    sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.data <- gene_summary
  }
  Seurat::VariableFeatures(sp.obj) <- hvgs
}
