#' Generate a plot of annotated spots.
#'
#' @name plotSpatialAttributes
#' @author Jack R. Leary
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object, plot the spatial location of each spot colored by a discrete attribute such as cluster ID, annotated domain, etc.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. Defaults to NULL.
#' @param attribute.plot A character specifying the categorical feature to be used to color the spots. Defaults to NULL.
#' @param pt.size A double specifying the size of the points to be plotted. Defaults to 2.
#' @param color.palette A vector containing colors that are passed to \code{\link[ggplot2]{scale_color_manual}} and thus used to color the spots. Defaults to NULL.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom dplyr select mutate
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom ggplot2 ggplot aes geom_point scale_y_continuous labs scale_color_manual
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link[Seurat]{SpatialDimPlot}}
#' @seealso \code{\link[ggspavis]{plotSpots}}
#' @seealso \code{\link{plotSpatialExpression}}
#' @export
#' @examples
#' data(seu_brain)
#' plotSpatialAttributes(seu_brain, attribute.plot = "region")

plotSpatialAttributes <- function(sp.obj = NULL,
                                  attribute.plot = NULL,
                                  pt.size = 2,
                                  color.palette = NULL) {
  # check inputs
  if (is.null(sp.obj) || is.null(attribute.plot)) { cli::cli_abort("Please provide all inputs to plotSpatialAttributes().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (inherits(sp.obj, "Seurat")) {
    meta_df <- sp.obj@meta.data
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    meta_df <- as.data.frame(sp.obj@metadata)
  }
  if (!attribute.plot %in% colnames(meta_df)) { cli::cli_abort("{.field attribute.plot} must exist in the colnames of the metadata of {.field sp.obj}.") }
  # generate plot
  if (inherits(sp.obj, "Seurat")) {
    coord_df <- Seurat::GetTissueCoordinates(sp.obj) %>%
                dplyr::select(1:2) %>%
                magrittr::set_colnames(c("x", "y")) %>%
                as.matrix() %>%
                scale() %>%
                as.data.frame() %>%
                dplyr::mutate(meta_vec = meta_df[, attribute.plot])
    p <- ggplot2::ggplot(coord_df, ggplot2::aes(x = y, y = x, color = meta_vec)) +
         ggplot2::geom_point(size = pt.size, stroke = 0) +
         ggplot2::scale_y_continuous(transform = "reverse")
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    coord_df <- SpatialExperiment::spatialCoords(sp.obj)
    coord1_name <- colnames(coord_df)[1]
    coord2_name <- colnames(coord_df)[2]
    p <- ggspavis::plotSpots(sp.obj,
                             x_coord = coord1_name,
                             y_coord = coord2_name,
                             annotate = attribute.plot,
                             point_size = pt.size,
                             show_axes = TRUE)
  }
  p <- p + ggplot2::labs(x = "Spatial 1",
                         y = "Spatial 2",
                         color = attribute.plot,
                         fill = attribute.plot)
  if (!is.null(color.palette)) {
    p <- p +
         ggplot2::scale_color_manual(values = color.palette)
  }
  p <- p +
       theme_bayesVG(spatial = TRUE)
  return(p)
}
