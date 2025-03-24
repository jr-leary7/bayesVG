#' Classify genes as spatially variable.
#'
#' @name classifySVGs
#' @author Jack R. Leary
#' @description After estimating per-gene spatial variation statistics, classify genes as spatially variable or not in a variety of ways.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. Defaults to NULL.
#' @param selection.method A string specifying what method should be used to classify genes as SVGs. Must be one of "rank", "quantile", or "cutoff". Defaults to "rank".
#' @param n.SVG An integer specifying the number of SVGs to select (if using rank-based selection). Defaults to 1000.
#' @param quantile.SVG A double specifying the quantile cutoff used to classify SVGs (if using quantile-based selection). Defaults to 0.75.
#' @param cutoff A double specifying the cutoff value for the spatial variation parameter \eqn{\tau_g} used to classify SVGs (if using cutoff-based selection). Defaults to 0.1.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom methods slot
#' @importFrom Seurat DefaultAssay VariableFeatures
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select arrange desc slice_head pull mutate filter if_else
#' @importFrom stats quantile
#' @importFrom S4Vectors DataFrame
#' @return An object of class \code{Seurat} with SVG metadata added.
#' @seealso \code{\link{findSpatiallyVariableFeaturesBayes}}
#' @seealso \code{\link[SeuratObject]{SVFInfo}}
#' @export
#' @examples
#' data(seu_brain)
#' seu_brain <- Seurat::SCTransform(seu_brain,
#'                                  assay = "Spatial",
#'                                  variable.features.n = 3000L,
#'                                  vst.flavor = "v2",
#'                                  return.only.var.genes = FALSE,
#'                                  seed.use = 312,
#'                                  verbose = FALSE)
#' seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain,
#'                                                 naive.hvgs = Seurat::VariableFeatures(seu_brain),
#'                                                 kernel = "matern",
#'                                                 kernel.smoothness = 1.5,
#'                                                 algorithm = "meanfield",
#'                                                 n.cores = 1L,
#'                                                 save.model = TRUE) %>%
#'              classifySVGs(n.SVG = 1000L)

classifySVGs <- function(sp.obj = NULL,
                         selection.method = "rank",
                         n.SVG = 1000L,
                         quantile.SVG = 0.75,
                         cutoff = 0.1) {
  # check inputs
  if (is.null(sp.obj)) { cli::cli_abort("Please provide an object to classifySVGs().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  selection.method <- tolower(selection.method)
  if (!selection.method %in% c("rank", "quantile", "cutoff")) { cli::cli_abort("Please provide a valid SVG selection method to classifySVGs().") }
  # extract gene spatial variation statistics
  if (inherits(sp.obj, "Seurat")) {
    version_check <- try({
      methods::slot(sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]], name = "meta.data")
    }, silent = TRUE)
    if (inherits(version_check, "try-error")) {
      gene_summary <- sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.features
    } else {
      gene_summary <- sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.data
    }
  } else {
    gene_summary <- SummarizedExperiment::rowData(sp.obj) %>%
                    as.data.frame()
  }
  # identify HVGs based on user-specified method
  if (selection.method == "rank") {
    svgs <- dplyr::arrange(gene_summary, amplitude_mean_rank) %>%
            dplyr::slice_head(n = n.SVG) %>%
            dplyr::pull(gene)
  } else if (selection.method == "quantile") {
    quantile_cutoff <- as.numeric(stats::quantile(gene_summary[, "amplitude_mean"], quantile.HVG))
    svgs <- dplyr::arrange(gene_summary, amplitude_mean_rank) %>%
            dplyr::filter(mean >= quantile_cutoff) %>%
            dplyr::pull(gene)
  } else if (selection.method == "cutoff") {
    svgs <- dplyr::arrange(gene_summary, amplitude_mean_rank) %>%
            dplyr::filter(amplitude_mean >= cutoff) %>%
            dplyr::pull(gene)
  }
  # add HVG classification back to object metadata
  gene_summary <- dplyr::mutate(gene_summary, svg_status = dplyr::if_else(gene %in% svgs, TRUE, FALSE))
  if (inherits(sp.obj, "Seurat")) {
    if (inherits(version_check, "try-error")) {
      sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.features <- gene_summary
    } else {
      sp.obj@assays[[Seurat::DefaultAssay(sp.obj)]]@meta.data <- gene_summary
    }
    Seurat::VariableFeatures(sp.obj) <- svgs
  } else {
    SummarizedExperiment::rowData(sp.obj) <- S4Vectors::DataFrame(gene_summary)
    sp.obj@metadata$svg_list <- svgs
  }
  return(sp.obj)
}
