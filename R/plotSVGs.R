#' Plot the results of SVG identification.
#'
#' @name plotSVGs
#' @author Jack R. Leary
#' @description This function plots the (log) mean versus the (log) dispersion of gene expression in order to visualize the relationship between the two statistics for HVGs and non-HVGs.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} upon which \code{\link{findSpatiallyVariableFeaturesBayes}} and \code{\link{classifySVGs}} have been run. Defaults to NULL.
#' @param pt.size A double specifying the size of the points on the plot. Defaults to 1.
#' @param pt.alpha A double specifying the opacity of the points on the plot. Defaults to 0.6.
#' @param add.smooth A Boolean specifying whether a GAM smoother should be overlaid in order to show the overall relationship between the (log) mean and (log) dispersion of gene expression. Defaults to TRUE.
#' @param n.genes.label An integer specifying the number of top HVGs to label on the plot using \code{\link[ggrepel]{geom_label_repel}}. If equal to 0, no genes will be labelled. Defaults to 10.
#' @param label.text.size A double specifying the size of the text for each top HVG's label. Defaults to 3.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom methods slot
#' @importFrom SummarizedExperiment rowData
#' @importFrom Seurat DefaultAssay
#' @importFrom dplyr inner_join mutate if_else arrange slice_head
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_color_manual labs guides guide_legend
#' @importFrom ggrepel geom_label_repel
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link{findVariableFeaturesBayes}}
#' @seealso \code{\link{classifyHVGs}}
#' @seealso \code{\link[Seurat]{VariableFeaturePlot}}
#' @export
#' @examples
#' data(seu_brain)
#' seu_brain <- Seurat::NormalizeData(seu_brain, verbose = FALSE) %>% 
#'              Seurat::FindVariableFeatures(nfeatures = 1000L, verbose = FALSE)
#' seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain,
#'                                                 naive.hvgs = Seurat::VariableFeatures(seu_brain),
#'                                                 kernel = "matern",
#'                                                 kernel.smoothness = 1.5,
#'                                                 algorithm = "meanfield",
#'                                                 n.cores = 1L,
#'                                                 save.model = TRUE) %>% 
#'              classifySVGs(n.SVG = 200L) 
#' plotSVGs(seu_brain)

plotSVGs <- function(sp.obj = NULL,
                     pt.size = 1,
                     pt.alpha = 0.6,
                     add.smooth = TRUE,
                     n.genes.label = 10L,
                     label.text.size = 3) {
  # check inputs
  if (is.null(sp.obj)) { cli::cli_abort("Please provide a {.pkg Seurat} or {.pkg SpatialExperiment} object to plotSVGs().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Argument {.field sp.obj} must be of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  # extract gene-level summary data.frame
  gene_summary <- getBayesianGeneStats(sp.obj, sort.values = FALSE)
  # generate naive gene statistics
  naive_stats <- computeNaiveGeneStatistics(sp.obj, use.norm = TRUE)
  genes_label <- dplyr::inner_join(gene_summary,
                                   naive_stats,
                                   by = "gene") %>%
                 dplyr::arrange(amplitude_mean_rank) %>%
                 dplyr::slice_head(n = n.genes.label)
  # generate plot
  p <- dplyr::inner_join(gene_summary,
                         naive_stats,
                         by = "gene") %>%
       dplyr::mutate(svg_label = dplyr::if_else(svg_status, "SVG", "Non-SVG")) %>%
       ggplot2::ggplot(ggplot2::aes(x = log(mu_naive), y = log(amplitude_mean), color = svg_label)) +
       ggplot2::geom_point(size = pt.size,
                           stroke = 0,
                           alpha = pt.alpha)
  if (add.smooth) {
    p <- p +
         ggplot2::geom_smooth(se = FALSE,
                              color = "dodgerblue",
                              method = "gam")
  }
  if (n.genes.label > 0) {
    p <- p +
         ggrepel::geom_label_repel(data = genes_label,
                                   mapping = ggplot2::aes(x = log(mu_naive), y = log(amplitude_mean), label = gene),
                                   inherit.aes = FALSE,
                                   force = 3,
                                   size = label.text.size,
                                   fontface = "italic")
  }
  p <- p +
       ggplot2::scale_color_manual(values = c("grey30", "firebrick")) +
       ggplot2::labs(x = expression(log(hat(mu)[g])),
                     y = expression(log(hat(tau)[g])),
                     color = "Status") +
       theme_bayesVG() +
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, stroke = 0.5, alpha = 1)))
  return(p)
}
