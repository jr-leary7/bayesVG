#' Generate a plot of annotated spots.
#'
#' @name plotSpatialExpression
#' @author Jack R. Leary
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object, plot the spatial location of each spot colored by the expression of a gene of interest.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. Defaults to NULL.
#' @param gene.plot A character specifying the gene whose expression will be used to color the spots. Defaults to NULL.
#' @param use.norm A Boolean specifying whether the raw or normalized counts should be plotted. Defaults to TRUE.
#' @param pt.size A double specifying the size of the points to be plotted. Defaults to 2.
#' @param color.palette A vector containing colors that are passed to \code{\link[ggplot2]{scale_color_gradientn}} and \code{\link[ggplot2]{scale_fill_gradientn}} and thus used to color the spots. Defaults to NULL.
#' @importFrom cli cli_abort
#' @importFrom Seurat GetAssayData DefaultAssay GetTissueCoordinates
#' @importFrom dplyr select mutate
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom ggplot2 ggplot aes geom_point scale_y_continuous labs scale_color_gradientn scale_fill_gradientn theme element_text
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link[Seurat]{SpatialFeaturePlot}}
#' @seealso \code{\link[ggspavis]{plotSpots}}
#' @export
#' @examples
#' data(seu_brain)
#' plotSpatialExpression(seu_brain,
#'                       gene.plot = "Nrgn",
#'                       use.norm = FALSE)



plotSpatialExpression <- function(sp.obj = NULL,
                                  gene.plot = NULL,
                                  use.norm = TRUE,
                                  pt.size = 2,
                                  color.palette = NULL) {
  # check inputs
  if (is.null(sp.obj) || is.null(gene.plot)) { cli::cli_abort("Please provide all inputs to plotSpatialExpression().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (!gene.plot %in% rownames(sp.obj)) { cli::cli_abort("{.field gene.plot} must exist in the rownames of {.field sp.obj}.") }
  # generate plot
  if (inherits(sp.obj, "Seurat")) {
    expr_vector <- Seurat::GetAssayData(sp.obj, assay = Seurat::DefaultAssay(sp.obj), layer = ifelse(use.norm, "data", "counts"))[gene.plot, ]
    coord_df <- Seurat::GetTissueCoordinates(seu_brain) %>% 
                dplyr::select(1:2) %>% 
                as.matrix() %>% 
                scale() %>% 
                as.data.frame() %>% 
                dplyr::mutate(gene_expr = expr_vector)
    p <- ggplot2::ggplot(coord_df, ggplot2::aes(x = y, y = x, color = gene_expr)) + 
         ggplot2::geom_point(size = pt.size, stroke = 0) + 
         ggplot2::scale_y_continuous(transform = "reverse")
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    if (!"gene" %in% colnames(SummarizedExperiment::rowData(sp.obj))) {
      SummarizedExperiment::rowData(sp.obj) <- as.data.frame(SummarizedExperiment::rowData(sp.obj)) %>%
                                               dplyr::mutate(gene = rownames(.), .before = 1) %>%
                                               S4Vectors::DataFrame()
    }
    coord_df <- SpatialExperiment::spatialCoords(sp.obj)
    coord1_name <- colnames(coord_df)[1]
    coord2_name <- colnames(coord_df)[2]
    p <- ggspavis::plotSpots(sp.obj,
                             x_coord = coord1_name,
                             y_coord = coord2_name,
                             annotate = gene.plot,
                             feature_names = "gene",
                             point_size = pt.size,
                             assay_name = ifelse(use.norm, "logcounts", "counts"),
                             show_axes = TRUE)
  }
  p <- p + ggplot2::labs(x = "Spatial 1", 
                         y = "Spatial 2", 
                         color = gene.plot, 
                         fill = gene.plot)
  if (!is.null(color.palette)) {
    p <- p + 
         ggplot2::scale_color_gradientn(colors = color.palette) + 
         ggplot2::scale_fill_gradientn(colors = color.palette)
  }
  p <- p + 
       theme_bayesVG(spatial = TRUE) + 
       ggplot2::theme(legend.title = ggplot2::element_text(face = "italic"))
  return(p)
}
