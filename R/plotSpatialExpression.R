#' Generate a plot of annotated spots.
#'
#' @name plotSpatialExpression
#' @author Jack R. Leary
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object, plot the spatial location of each spot colored by the expression of a gene of interest.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. Defaults to NULL.
#' @param gene.plot A character or character vector specifying the gene(s) whose expression will be used to color the spots. Defaults to NULL.
#' @param use.norm A Boolean specifying whether the raw or normalized counts should be plotted. Defaults to TRUE.
#' @param scale.expression A Boolean specifying whether expression should be scaled prior to plotting. Mostly relevant for when \code{gene.plot} contains more than one value. Defaults to TRUE when \code{length(gene.plot) > 1} and FALSE otherwise.
#' @param pt.size A double specifying the size of the points to be plotted. Defaults to 2.
#' @param color.palette A vector containing colors that are passed to \code{\link[ggplot2]{scale_color_gradientn}} and thus used to color the spots. Defaults to NULL.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom Seurat GetAssayData DefaultAssay GetTissueCoordinates
#' @importFrom Matrix t
#' @importFrom tidyr pivot_longer
#' @importFrom coop scaler
#' @importFrom dplyr select bind_cols mutate
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom BiocGenerics counts
#' @importFrom ggplot2 ggplot aes facet_wrap geom_point scale_y_continuous labs scale_color_gradientn theme element_text element_blank element_rect
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link[Seurat]{SpatialFeaturePlot}}
#' @export
#' @examples
#' data(seu_brain)
#' plotSpatialExpression(seu_brain,
#'                       gene.plot = "Nrgn",
#'                       use.norm = FALSE)



plotSpatialExpression <- function(sp.obj = NULL,
                                  gene.plot = NULL,
                                  use.norm = TRUE,
                                  scale.expression = NULL, 
                                  pt.size = 2,
                                  color.palette = NULL) {
  # check inputs
  if (is.null(sp.obj) || is.null(gene.plot)) { cli::cli_abort("Please provide all inputs to plotSpatialExpression().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (!all(gene.plot %in% rownames(sp.obj))) { cli::cli_abort("All values of {.field gene.plot} must exist in the rownames of {.field sp.obj}.") }
  if (is.null(scale.expression)) {
    if (length(gene.plot) > 1) {
      scale.expression <- TRUE
    } else {
      scale.expression <- FALSE
    }
  }
  # extract spatial coordinates & gene expression from sp.obj
  if (inherits(sp.obj, "Seurat")) {
    expr_mtx <- Seurat::GetAssayData(sp.obj,
                                     assay = Seurat::DefaultAssay(sp.obj),
                                     layer = ifelse(use.norm, "data", "counts"))[gene.plot, , drop = FALSE]
    expr_df <- as.data.frame(Matrix::t(expr_mtx))
    coord_df <- Seurat::GetTissueCoordinates(sp.obj) %>%
                dplyr::select(1:2) %>%
                as.matrix() %>%
                coop::scaler() %>%
                as.data.frame() %>%
                magrittr::set_colnames(c("x", "y")) %>%
                dplyr::bind_cols(expr_df)
    if (length(gene.plot) > 1) {
      coord_df <- tidyr::pivot_longer(coord_df, 
                                      cols = !c(x, y), 
                                      names_to = "gene", 
                                      values_to = "gene_expr")
    } else {
      colnames(coord_df)[length(colnames(coord_df))] <- "gene_expr"
    }
    if (scale.expression) {
      coord_df <- dplyr::mutate(coord_df, gene_expr = as.numeric(coop::scaler(gene_expr)))
    }
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    if (use.norm) {
      expr_mtx <- SingleCellExperiment::logcounts(sp.obj)[gene.plot, , drop = FALSE]
    } else {
      expr_mtx <- BiocGenerics::counts(sp.obj)[gene.plot, , drop = FALSE]
    }
    expr_df <- as.data.frame(Matrix::t(expr_mtx))
    coord_df <- SpatialExperiment::spatialCoords(sp.obj) %>% 
                coop::scaler() %>% 
                as.data.frame() %>% 
                magrittr::set_colnames(c("x", "y")) %>%
                dplyr::bind_cols(expr_df)
    if (length(gene.plot) > 1) {
      coord_df <- tidyr::pivot_longer(coord_df, 
                                      cols = !c(x, y), 
                                      names_to = "gene", 
                                      values_to = "gene_expr")
    } else {
      colnames(coord_df)[length(colnames(coord_df))] <- "gene_expr"
    }
    if (scale.expression) {
      coord_df <- dplyr::mutate(coord_df, gene_expr = as.numeric(coop::scaler(gene_expr)))
    }
  }
  # generate (potentially facetted) plot along with labels and colors
  if (length(gene.plot) > 1) {
    p <- ggplot2::ggplot(coord_df, ggplot2::aes(x = y, y = x, color = gene_expr)) +
         ggplot2::facet_wrap(~gene, scales = "fixed") + 
         ggplot2::geom_point(size = pt.size, stroke = 0) +
         ggplot2::scale_y_continuous(transform = "reverse")
  } else {
    p <- ggplot2::ggplot(coord_df, ggplot2::aes(x = y, y = x, color = gene_expr)) +
         ggplot2::geom_point(size = pt.size, stroke = 0) +
         ggplot2::scale_y_continuous(transform = "reverse")
  }
  p <- p + ggplot2::labs(x = "Spatial 1",
                         y = "Spatial 2",
                         color = ifelse(length(gene.plot) > 1, 
                                        "Scaled\nExpression", 
                                        gene.plot))
  if (!is.null(color.palette)) {
    p <- p +
         ggplot2::scale_color_gradientn(colors = color.palette)
  }
  p <- p +
       theme_bayesVG(spatial = TRUE) +
       ggplot2::theme(plot.title = ggplot2::element_blank())
  if (length(gene.plot) > 1) {
    p <- p + 
         ggplot2::theme(strip.text = ggplot2::element_text(face = "italic"), 
                        strip.clip = "on",
                        strip.background = ggplot2::element_rect(linewidth = 2 * theme_bayesVG()$line$linewidth))
  } else {
    p <- p + 
         ggplot2::theme(legend.title = ggplot2::element_text(face = "italic"))
  }
  return(p)
}
