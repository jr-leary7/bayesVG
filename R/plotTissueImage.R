#' Generate a plot of the tissue image given a \code{Seurat} object.
#'
#' @name plotTissueImage
#' @author Jack R. Leary
#' @description Given a \code{Seurat} object, plot a high-definition visual of the tissue image without any annotations.
#' @param seu.obj An object of class \code{Seurat} containing spatial data. Defaults to NULL.
#' @param image.name A string specifying the name of the image to extract from \code{seu.obj}. If left unspecified, an attempt will be made to select the correct image automatically. Defaults to NULL. 
#' @import magrittr 
#' @importFrom cli cli_abort
#' @importFrom Seurat GetImage
#' @importFrom dplyr mutate
#' @importFrom grDevices rgb
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_identity scale_y_continuous coord_fixed labs
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link[Seurat]{GetImage}}
#' @seealso \code{\link[Seurat]{SpatialFeaturePlot}}
#' @export
#' @examples
#' data(seu_brain)
#' plotTissueImage(seu_brain)

plotTissueImage <- function(seu.obj = NULL, image.name = NULL) {
  # check inputs 
  if (is.null(seu.obj)) { cli::cli_abort("Please provide all inputs to plotTissueImage().") }
  # extract image array & reorder 
  img <- Seurat::GetImage(seu.obj, 
                          mode = "raw", 
                          image = image.name)
  img_xy <- aperm(img, c(2, 1, 3)) 
  # generate tissue image plot 
  img_df <- expand.grid(x = seq_len(dim(img)[2]), y = seq_len(dim(img)[1])) %>%
            dplyr::mutate(r = as.vector(img_xy[, , 1]),
                          g = as.vector(img_xy[, , 2]),
                          b = as.vector(img_xy[, , 3]),
                          hex = grDevices::rgb(r, g, b))
  p <- ggplot2::ggplot(img_df, ggplot2::aes(x = x, y = y, fill = hex)) +
       ggplot2::geom_raster() +
       ggplot2::scale_fill_identity() + 
       ggplot2::scale_y_continuous(trans = "reverse") + 
       ggplot2::coord_fixed(expand = FALSE) + 
       ggplot2::labs(x = "Spatial 1", y = "Spatial 2") + 
       theme_bayesVG(spatial = TRUE)
  return(p)
}
