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
#' @import magrittr
#' @importFrom methods slot
#' @importFrom SummarizedExperiment rowData
#' @importFrom Seurat DefaultAssay
#' @importFrom dplyr mutate if_else arrange desc slice_head
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth scale_color_manual labs guides guide_legend
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
#'             classifyHVGs(n.HVG = 1000L)
#' plotHVGs(seu_pbmc)

plotHVGs <- function(sc.obj = NULL,
                     pt.size = 1,
                     pt.alpha = 0.6,
                     add.smooth = TRUE,
                     n.genes.label = 10L,
                     label.text.size = 3) {
  # check inputs
  if (is.null(sc.obj)) { stop("Please provide a Seurat or SingleCellExperiment object to plotHVGs().") }
  if (!(inherits(sc.obj, "Seurat") || inherits(sc.obj, "SingleCellExperiment"))) { stop("Argument sc.obj must be of class Seurat or SingleCellExperiment.") }
  # extract gene-level summary data.frame
  if (inherits(sc.obj, "SingleCellExperiment")) {
    gene_summary <- as.data.frame(SummarizedExperiment::rowData(sc.obj))
  } else if (inherits(sc.obj, "Seurat")) {
    version_check <- try({
      methods::slot(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]], name = "meta.data")
    }, silent = TRUE)
    if (inherits(version_check, "try-error")) {
      gene_summary <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features
    } else {
      gene_summary <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data
    }
  }
  genes_label <- dplyr::arrange(gene_summary, dplyr::desc(dispersion_mean)) %>%
                 dplyr::slice_head(n = n.genes.label)
  # generate plot
  p <- dplyr::mutate(gene_summary, hvg_label = dplyr::if_else(hvg_status, "HVG", "Non-HVG")) %>%
       ggplot2::ggplot(ggplot2::aes(x = log(mu_mean), y = log(dispersion_mean), color = hvg_label)) +
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
                                   mapping = aes(x = log(mu_mean), y = log(dispersion_mean), label = gene),
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
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, stroke = 0.5, alpha = 1)))
  return(p)
}
