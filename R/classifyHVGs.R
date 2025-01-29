#' Classify genes as highly variable.
#'
#' @name classifyHVGs
#' @author Jack R. Leary
#' @description After estimating per-gene statistics, classify genes as highly variable or not in a variety of ways.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @param selection.variable A string specifying the random variable used to select HVGs. Must be one of "dispersion" or "sigma2". Defaults to "dispersion".
#' @param selection.method A string specifying what method should be used to classify genes as HVGs. Must be one of "rank", "quantile", or "cutoff". Defaults to "rank".
#' @param n.HVG An integer specifying the number of HVGs to select (if using rank-based selection). Defaults to 2000.
#' @param quantile.HVG A double specifying the quantile cutoff used to classify HVGs (if using quantile-based selection). Defaults to 0.75.
#' @param cutoff A double specifying the cutoff value for dispersion or variance (depending on how \code{selection.variable} is defined) used to classify HVGs (if using cutoff-based selection). Defaults to 3.
#' @import magrittr
#' @importFrom rlang sym
#' @importFrom SingleCellExperiment rowData
#' @importFrom Seurat DefaultAssay VariableFeatures
#' @importFrom dplyr select arrange desc slice_head pull mutate filter if_else
#' @importFrom stats quantile
#' @importFrom S4Vectors DataFrame
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added.
#' @seealso \code{\link{findVariableFeaturesBayes}}
#' @seealso \code{\link[SeuratObject]{HVFInfo}}
#' @export

classifyHVGs <- function(sc.obj = NULL,
                         selection.variable = "dispersion",
                         selection.method = "rank",
                         n.HVG = 2000L,
                         quantile.HVG = 0.75,
                         cutoff = 3) {
  # check inputs
  if (is.null(sc.obj)) { stop("Please provide an object to classifyHVGs().") }
  selection.variable <- tolower(selection.variable)
  if (!selection.variable %in% c("dispersion", "sigma2")) { stop("Please provide a valid random variable used in classifyHVGs().") }
  selection.method <- tolower(selection.method)
  if (!selection.method %in% c("rank", "quantile", "cutoff")) { stop("Please provide a valid HVG selection method to classifyHVGs().") }
  # set up ranking variable
  ranker_var <- rlang::sym(paste0(selection.variable, "_mean"))
  # extract gene mean & dispersion statistics
  if (inherits(sc.obj, "SingleCellExperiment")) {
    gene_summary <- as.data.frame(SingleCellExperiment::rowData(sc.obj))
  } else if (inherits(sc.obj, "Seurat")) {
    version_check <- try({
      slot(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]], name = "meta.data")
    }, silent = TRUE)
    if (inherits(version_check, "try-error")) {
      gene_summary <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features
    } else {
      gene_summary <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data
    }
  }
  # identify HVGs based on user-specified method
  if (selection.method == "rank") {
    hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
            dplyr::slice_head(n = n.HVG) %>%
            dplyr::pull(gene)
  } else if (selection.method == "quantile") {
    quantile_cutoff <- as.numeric(stats::quantile(gene_summary[, paste0(selection.variable, "_mean")], quantile.HVG))
    hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
            dplyr::filter(!!ranker_var >= quantile_cutoff) %>%
            dplyr::pull(gene)
  } else if (selection.method == "cutoff") {
    hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
            dplyr::filter(!!ranker_var >= cutoff) %>%
            dplyr::pull(gene)
  }
  # add HVG classification back to object metadata
  gene_summary <- dplyr::mutate(gene_summary, hvg = dplyr::if_else(gene %in% hvgs, TRUE, FALSE))
  if (inherits(sc.obj, "SingleCellExperiment")) {
    SingleCellExperiment::rowData(sc.obj) <- S4Vectors::DataFrame(gene_summary)
  } else if (inherits(sc.obj, "Seurat")) {
    if (inherits(version_check, "try-error")) {
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features <- gene_summary
    } else {
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- gene_summary
    }
    Seurat::VariableFeatures(sc.obj) <- hvgs
  }
  return(sc.obj)
}
