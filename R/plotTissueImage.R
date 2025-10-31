#' Generate a plot of the tissue image without any annotations.
#'
#' @name plotTissueImage
#' @author Jack R. Leary
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object, plot a high-definition visual of the tissue image without any annotations.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial image data. Defaults to NULL.
#' @param image.name A string specifying the name of the image to extract when \code{sp.obj} is of class \code{Seurat}. If left unspecified, an attempt will be made to select the correct image automatically. Defaults to NULL. 
#' @import magrittr 
#' @importFrom cli cli_abort
#' @importFrom Seurat GetImage
#' @importFrom dplyr mutate
#' @importFrom grDevices rgb
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_identity scale_y_continuous coord_fixed labs
#' @importFrom SpatialExperiment imgRaster
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link[Seurat]{GetImage}}
#' @seealso \code{\link[Seurat]{SpatialFeaturePlot}}
#' @seealso \code{\link[SpatialExperiment]{imgRaster}}
#' @export
#' @examples
#' data(seu_brain)
#' plotTissueImage(seu_brain)

plotTissueImage <- function(sp.obj = NULL, image.name = NULL) {
  # check inputs 
  if (is.null(sp.obj)) { cli::cli_abort("Please provide all inputs to plotTissueImage().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (inherits(sp.obj, "Seurat")) {
    if (!is.null(image.name) && !image.name %in% names(sp.obj@images)) { 
      cli::cli_abort("Please provide a valid image name.") 
    }
  }
  if (inherits(sp.obj, "Seurat")) {
    # extract image array & reorder 
    img <- Seurat::GetImage(sp.obj, 
                            mode = "raw", 
                            image = image.name)
    img_xy <- aperm(img, perm = c(2, 1, 3)) 
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
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    # extract image data & convert to data.frame
    img <- SpatialExperiment::imgRaster(sp.obj)
    img_mat <- as.matrix(img)
    img_mat <- substr(img_mat, 
                      start = 1, 
                      stop = 7)
    img_df <- expand.grid(x = seq_len(ncol(img_mat)), y = seq_len(nrow(img_mat)))
    img_df$fill <- as.vector(img_mat)
    # generate tissue image plot 
    p <- ggplot2::ggplot(img_df, aes(x = x, y = y, fill = fill)) + 
         ggplot2::geom_raster() +
         ggplot2::scale_fill_identity() + 
         ggplot2::scale_y_continuous(trans = "reverse") + 
         ggplot2::coord_fixed(expand = FALSE) + 
         ggplot2::labs(x = "Spatial 1", y = "Spatial 2") + 
         theme_bayesVG(spatial = TRUE)
  }
  return(p)
}
