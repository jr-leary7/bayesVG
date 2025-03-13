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
#' @importFrom cli cli_abort
#' @importFrom rlang sym
#' @importFrom purrr reduce
#' @importFrom SummarizedExperiment rowData
#' @importFrom Seurat DefaultAssay VariableFeatures
#' @importFrom dplyr select arrange desc slice_head pull mutate filter if_else
#' @importFrom stats quantile
#' @importFrom methods slot
#' @importFrom S4Vectors DataFrame
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with HVG metadata added.
#' @seealso \code{\link{findVariableFeaturesBayes}}
#' @seealso \code{\link[SeuratObject]{HVFInfo}}
#' @export
#' @examples
#' data(seu_pbmc)
#' seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
#'                                       n.cells.subsample = 500L,
#'                                       algorithm = "meanfield",
#'                                       n.cores.per.chain = 1L,
#'                                       save.model = TRUE) %>%
#'             classifyHVGs(n.HVG = 3000L)

classifyHVGs <- function(sc.obj = NULL,
                         selection.variable = "dispersion",
                         selection.method = "rank",
                         n.HVG = 2000L,
                         quantile.HVG = 0.75,
                         cutoff = 3) {
  # check inputs
  if (is.null(sc.obj)) { cli::cli_abort("Please provide an object to classifyHVGs().") }
  selection.variable <- tolower(selection.variable)
  if (!selection.variable %in% c("dispersion", "sigma2")) { cli::cli_abort("Please provide a valid random variable used in classifyHVGs().") }
  selection.method <- tolower(selection.method)
  if (!selection.method %in% c("rank", "quantile", "cutoff")) { cli::cli_abort("Please provide a valid HVG selection method to classifyHVGs().") }
  # set up ranking variable
  ranker_var <- rlang::sym(paste0(selection.variable, "_mean"))
  # check if data are single- or multi-subject
  if (inherits(sc.obj, "SingleCellExperiment")) {
    multi_subject_flag <- ifelse(is.null(sc.obj@metadata$gene_stats_bayes), 
                                 FALSE, 
                                 TRUE)
  } else if (inherits(sc.obj, "Seurat")) {
    multi_subject_flag <- ifelse(is.null(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@misc$gene_stats_bayes), 
                                 FALSE, 
                                 TRUE)
  }
  # extract gene mean & dispersion statistics
  if (inherits(sc.obj, "SingleCellExperiment")) {
    if (multi_subject_flag) {
      gene_summary <- purrr::reduce(sc.obj@metadata$gene_stats_bayes, rbind)
    } else {
      gene_summary <- as.data.frame(SummarizedExperiment::rowData(sc.obj))
    }
  } else if (inherits(sc.obj, "Seurat")) {
    if (multi_subject_flag) {
      gene_summary <- purrr::reduce(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@misc$gene_stats_bayes, rbind)
    } else {
      version_check <- try({
        methods::slot(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]], name = "meta.data")
      }, silent = TRUE)
      if (inherits(version_check, "try-error")) {
        gene_summary <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features
      } else {
        gene_summary <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data
      }
    }
  }
  # identify HVGs based on user-specified method
  if (selection.method == "rank") {
    if (multi_subject_flag) {
      hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
              dplyr::distinct(gene) %>%
              dplyr::slice_head(n = n.HVG) %>%
              dplyr::pull(gene)
    } else {
      hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
              dplyr::slice_head(n = n.HVG) %>%
              dplyr::distinct(gene) %>%
              dplyr::pull(gene)
    }
  } else if (selection.method == "quantile") {
    quantile_cutoff <- as.numeric(stats::quantile(gene_summary[, paste0(selection.variable, "_mean")], quantile.HVG))
    hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
            dplyr::filter(!!ranker_var >= quantile_cutoff) %>%
            dplyr::distinct(gene) %>%
            dplyr::pull(gene)
  } else if (selection.method == "cutoff") {
    hvgs <- dplyr::arrange(gene_summary, dplyr::desc(!!ranker_var)) %>%
            dplyr::filter(!!ranker_var >= cutoff) %>%
            dplyr::distinct(gene) %>%
            dplyr::pull(gene)
  }
  # add HVG classification back to object metadata
  gene_summary <- dplyr::mutate(gene_summary, hvg_status = dplyr::if_else(gene %in% hvgs, TRUE, FALSE))
  if (inherits(sc.obj, "SingleCellExperiment")) {
    if (multi_subject_flag) {
      gene_summary_list <- split(gene_summary, gene_summary$subject)
      gene_summary_list <- purrr::map(gene_summary_list, \(x) {
        rownames(x) <- x$gene
        x <- x[rownames(sc.obj), ]
        return(x)
      })
      sc.obj@metadata$gene_stats_bayes <- gene_summary_list
    } else {
      SummarizedExperiment::rowData(sc.obj) <- S4Vectors::DataFrame(gene_summary)
    }
  } else if (inherits(sc.obj, "Seurat")) {
    if (multi_subject_flag) {
      gene_summary_list <- split(gene_summary, gene_summary$subject)
      gene_summary_list <- purrr::map(gene_summary_list, \(x) {
        rownames(x) <- x$gene
        x <- x[rownames(sc.obj), ]
        return(x)
      })
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@misc$gene_stats_bayes <- gene_summary_list
    } else {
      version_check <- try({
        methods::slot(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]], name = "meta.data")
      }, silent = TRUE)
      if (inherits(version_check, "try-error")) {
        sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features <- gene_summary
      } else {
        sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- gene_summary
      }
    }
    Seurat::VariableFeatures(sc.obj) <- hvgs
  }
  return(sc.obj)
}
