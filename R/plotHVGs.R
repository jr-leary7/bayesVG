#' Plot the results of HVG identification.
#'
#' @name plotHVGs
#' @author Jack R. Leary
#' @description This function plots the (log) mean versus the (log) dispersion of gene expression in order to visualize the relationship between the two statistics for HVGs and non-HVGs.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment} upon which \code{\link{findVariableFeaturesBayes}} and \code{\link{classifyHVGs}} have been run. Defaults to NULL.
#' @param pt.size A double specifying the size of the points on the plot. Defaults to 1.
#' @param pt.alpha A double specifying the opacity of the points on the plot. Defaults to 0.6.
#' @param add.smooth A Boolean specifying whether a GAM smoother should be overlaid in order to show the overall relationship between the (log) mean and (log) dispersion of gene expression. Defaults to TRUE.
#' @param n.genes.label An integer specifying the number of top HVGs to label on the plot using \code{\link[ggrepel]{geom_label_repel}}. If equal to 0, no genes will be labelled. Defaults to 10.
#' @param label.text.size A double specifying the size of the text for each top HVG's label. Defaults to 3.
#' @details
#' \itemize{
#' \item When multiple subjects are present (this is detected automatically by internal logic), the resulting plot will be faceted by subject and the top \code{n.genes.label} per-subject will be labelled. 
#' }
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom methods slot
#' @importFrom SummarizedExperiment rowData
#' @importFrom Seurat DefaultAssay
#' @importFrom dplyr mutate if_else arrange desc slice_head with_groups
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_smooth scale_color_manual labs guides guide_legend
#' @importFrom ggrepel geom_label_repel
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link{findVariableFeaturesBayes}}
#' @seealso \code{\link{classifyHVGs}}
#' @seealso \code{\link[Seurat]{VariableFeaturePlot}}
#' @export
#' @examples
#' data(seu_pbmc)
#' seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
#'                                       n.cells.subsample = 500L,
#'                                       algorithm = "meanfield",
#'                                       n.cores.per.chain = 1L,
#'                                       save.model = TRUE) %>%
#'             classifyHVGs(n.HVG = 500L)
#' plotHVGs(seu_pbmc)

plotHVGs <- function(sc.obj = NULL,
                     pt.size = 1,
                     pt.alpha = 0.6,
                     add.smooth = TRUE,
                     n.genes.label = 10L,
                     label.text.size = 3) {
  # check inputs
  if (is.null(sc.obj)) { cli::cli_abort("Please provide a {.pkg Seurat} or {.pkg SingleCellExperiment} object to plotHVGs().") }
  if (!(inherits(sc.obj, "Seurat") || inherits(sc.obj, "SingleCellExperiment"))) { cli::cli_abort("Argument {.field sc.obj} must be of class {.pkg Seurat} or {.pkg SingleCellExperiment}.") }
  # extract gene-level summary data.frame
  gene_summary <- getBayesianGeneStats(sc.obj, sort.values = FALSE)
  multi_subject_flag <- ifelse("subject" %in% colnames(gene_summary), 
                               TRUE, 
                               FALSE)
  if (multi_subject_flag) {
    genes_label <- dplyr::arrange(gene_summary, dplyr::desc(dispersion_mean)) %>%
                   dplyr::with_groups(subject, 
                                      dplyr::slice_head, 
                                      n = 10)
  } else {
    genes_label <- dplyr::arrange(gene_summary, dplyr::desc(dispersion_mean)) %>%
                   dplyr::slice_head(n = n.genes.label)
  }
  # generate plot
  p <- dplyr::mutate(gene_summary, hvg_label = dplyr::if_else(hvg_status, "HVG", "Non-HVG")) %>%
       ggplot2::ggplot(ggplot2::aes(x = log(mu_mean), y = log(dispersion_mean), color = hvg_label)) +
       ggplot2::geom_point(size = pt.size,
                           stroke = 0,
                           alpha = pt.alpha)
  if (multi_subject_flag) {
    p <- p + ggplot2::facet_wrap(subject)
  }
  if (add.smooth) {
    p <- p +
         ggplot2::geom_smooth(se = FALSE,
                              color = "dodgerblue",
                              method = "gam")
  }
  if (n.genes.label > 0) {
    p <- p +
         ggrepel::geom_label_repel(data = genes_label,
                                   mapping = ggplot2::aes(x = log(mu_mean), y = log(dispersion_mean), label = gene),
                                   inherit.aes = FALSE,
                                   force = 3,
                                   size = label.text.size,
                                   fontface = "italic")
}
  p <- p +
       ggplot2::scale_color_manual(values = c("firebrick", "grey30")) +
       ggplot2::labs(x = expression(log(hat(mu)[g])),
                     y = expression(log(hat(d)[g])),
                     color = "Status") +
       theme_bayesVG() +
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, 
                                                                         stroke = 0.5, 
                                                                         alpha = 1)))
  return(p)
}
