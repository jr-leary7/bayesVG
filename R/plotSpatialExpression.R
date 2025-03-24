#' Generate a plot of annotated spots. 
#' 
#' @name plotSpatialExpression
#' @author Jack R. Leary 
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object, plot the spatial location of each spot colored by the expression of a gene of interest.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. Defaults to NULL. 
#' @param gene.plot A character specifying the gene whose expression will be used to color the spots. Defaults to NULL. 
#' @param use.norm A Boolean specifying whether the raw or normalized counts should be plotted. Defaults to TRUE. 
#' @param color.palette A vector containing colors that are passed to \code{\link[ggplot2]{scale_color_gradientn}} and thus used to color the spots. Defaults to NULL.
#' @importFrom cli cli_abort
#' @importFrom Seurat SpatialFeaturePlot
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr mutate
#' @importFrom S4Vectors DataFrame
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom ggplot2 labs scale_color_gradientn
#' @return An object of class \code{ggplot2}. 
#' @seealso \code{\link[Seurat]{SpatialFeaturePlot}}
#' @seealso \code{\link[ggspavis]{plotSpots}}
#' @export 
#' @examples 
#' data(seu_brain)
#' plotSpatialExpression(seu_brain, gene.plot = "Nrgn")



plotSpatialExpression <- function(sp.obj = NULL, 
                                  gene.plot = NULL, 
                                  use.norm = TRUE, 
                                  color.palette = NULL) {
  # check inputs 
  if (is.null(sp.obj) || is.null(gene.plot)) { cli::cli_abort("Please provide all inputs to spotPlot().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class Seurat or SpatialExperiment.") }
  if (!gene.plot %in% rownames(sp.obj)) { cli::cli_abort("gene.plot must exist in the rownames of sp.obj.") }
  # generate plot 
  if (inherits(sp.obj, "Seurat")) {
    p <- Seurat::SpatialFeaturePlot(sp.obj, 
                                    features = gene.plot, 
                                    slot = ifelse(use.norm, "data", "counts"))
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    if (!"gene" %in% colnames(SummarizedExperiment::rowData(sp.obj))) {
      SummarizedExperiment::rowData(sp.obj) <- as.data.frame(SummarizedExperiment::rowData(sp.obj)) %>% 
                                               dplyr::mutate(gene = rownames(.), .before = 1) %>% 
                                               S4Vectors::DataFrame()
    }
    coord_df <- SpatialExperiment::spatialCoords(sp.obj)
    coord1_name <- colnames(coord_df)[2]
    coord2_name <- colnames(coord_df)[1]
    p <- ggspavis::plotSpots(sp.obj, 
                             x_coord = coord1_name, 
                             y_coord = coord2_name,
                             annotate = gene.plot, 
                             feature_names = "gene", 
                             point_size = 2, 
                             assay_name = ifelse(use.norm, "logcounts", "counts"), 
                             show_axes = TRUE)
  }
  p <- p + ggplot2::labs(x = "Spatial 1", y = "Spatial 2")
  if (!is.null(color.palette)) {
    p <- p + ggplot2::scale_color_gradientn(colors = color.palette)
  }
  p <- p + theme_bayesVG(spatial = TRUE)
  return(p)
}
