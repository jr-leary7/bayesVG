#' Generate a plot of annotated spots.
#'
#' @name plotModuleScores
#' @author Jack R. Leary
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object, generate a variety of plots showing the distribution of spatial module scores.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. Defaults to NULL.
#' @param module.plot A character specifying the name of the spatial module whose values should be visualized. Use e.g., "1" to plot module 1. Defaults to NULL.
#' @param plot.type A character specifying which type of plot should be generated. Must be one of "spatial", "embedding", or "violin". Defaults to "spatial".
#' @param embedding.name A character specifying which low-dimensional embedding to be used when \code{plot.type = "embedding"}. Defaults to NULL.
#' @param violin.group A character specifying which categorical variable in the metadata of \code{sp.obj} will be used to group the violins when \code{plot.type = "violin"}. Defaults to NULL.
#' @param pt.size A double specifying the size of the points to be plotted when \code{plot.type} is "spatial" or "embedding". Defaults to 2.
#' @param pt.alpha A double specifying the opacity of the points to be plotted when \code{plot.type} is "spatial" or "embedding". Defaults to 0.75.
#' @param color.palette A vector containing colors that are passed to \code{\link[ggplot2]{scale_color_gradientn}} or \code{\link[ggplot2]{scale_color_manual}} and \code{\link[ggplot2]{scale_fill_manual}}, depending on the value of \code{plot.type}. Defaults to NULL.
#' @import magrittr
#' @importFrom cli cli_abort
#' @importFrom rlang sym
#' @importFrom Seurat Reductions GetTissueCoordinates Embeddings
#' @importFrom dplyr select mutate
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom ggplot2 ggplot aes geom_point geom_violin scale_y_continuous labs scale_color_gradientn scale_color_manual scale_fill_manual
#' @importFrom utils packageVersion
#' @return An object of class \code{ggplot2}.
#' @seealso \code{\link{scoreSpatialModules}}
#' @seealso \code{\link{plotSpatialExpression}}
#' @seealso \code{\link{plotSpatialAttributes}}
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
#' svg_clusters <- clusterSVGsBayes(seu_brain,
#'                                  svgs = Seurat::VariableFeatures(seu_brain),
#'                                  n.clust = 2L,
#'                                  n.cores = 1L)
#' seu_brain <- scoreSpatialModules(seu_brain,
#'                                  svg.clusters = svg_clusters,
#'                                  n.cores = 1L)
#' plotModuleScores(seu_brain, module.plot = "1")



plotModuleScores <- function(sp.obj = NULL,
                             module.plot = NULL,
                             plot.type = "spatial",
                             embedding.name  = NULL,
                             violin.group = NULL,
                             pt.size = 2,
                             pt.alpha = 0.75,
                             color.palette = NULL) {
  # check inputs
  if (is.null(sp.obj) || is.null(module.plot)) { cli::cli_abort("Please provide all inputs to plotModuleScores().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (inherits(sp.obj, "Seurat")) {
    meta_df <- sp.obj@meta.data
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    meta_df <- as.data.frame(sp.obj@metadata)
  }
  module_name <- paste0("svg_cluster_", module.plot, "_UCell")
  module_name_sym <- rlang::sym(module_name)
  module_label <- paste0("Module ", module.plot, " Score")
  if (!module_name %in% colnames(meta_df)) { cli::cli_abort("Please provide a valid module name; check the metadata of {.field sp.obj}. to make sure the specified module exists.") }
  if (plot.type == "embedding") {
    if (inherits(sp.obj, "Seurat")) {
      embedding_names <- Seurat::Reductions(sp.obj)
    } else if (inherits(sp.obj, "SpatialExperiment")) {
      embedding_names <- SingleCellExperiment::reducedDimNames(sp.obj)
    }
    if (!embedding.name %in% embedding_names) { cli::cli_abort("{.field embedding.name} must be the name of a low-dimensional embedding present in {.field sp.obj}. Check to make sure you've spelled / capitalized {.field embedding.name} correctly.") }
  } else if (plot.type == "violin") {
    violin_group_sym <- rlang::sym(violin.group)
    if (!violin.group %in% colnames(meta_df)) { cli::cli_abort("{.field violin.group} must exist in the colnames of the metadata of {.field sp.obj}.") }
  }
  # generate plot
  if (plot.type == "spatial") {
    if (inherits(sp.obj, "Seurat")) {
      coord_df <- Seurat::GetTissueCoordinates(sp.obj)
    } else if (inherits(sp.obj, "SpatialExperiment")) {
      coord_df <- SpatialExperiment::spatialCoords(sp.obj)
    }
    coord_df <- dplyr::select(coord_df, 1:2) %>%
                magrittr::set_colnames(c("x", "y")) %>%
                as.matrix() %>%
                scale() %>%
                as.data.frame() %>%
                dplyr::mutate(module_score = meta_df[[module_name]])
    p <- ggplot2::ggplot(coord_df, ggplot2::aes(x = y, y = x, color = module_score)) +
         ggplot2::geom_point(size = pt.size,
                             alpha = pt.alpha,
                             stroke = 0) +
         ggplot2::scale_y_continuous(transform = "reverse") +
         ggplot2::labs(x = "Spatial 1",
                       y = "Spatial 2",
                       color = module_label) +
         theme_bayesVG(spatial = TRUE)
    if (!is.null(color.palette)) {
      p <- p + ggplot2::scale_color_gradientn(colours = color.palette)
    }
  } else if (plot.type == "embedding") {
    if (inherits(sp.obj, "Seurat")) {
      embedding_df <- Seurat::Embeddings(sp.obj, reduction = embedding.name)
    } else if (inherits(sp.obj, "SpatialExperiment")) {
      embedding_df <- SingleCellExperiment::reducedDim(sp.obj, type = embedding.name)
    }
    embedding_labels <- paste0(toupper(embedding.name), " ", 1:2)
    embedding_df <- dplyr::select(as.data.frame(embedding_df), 1:2) %>%
                    magrittr::set_colnames(c("dim1", "dim2")) %>%
                    dplyr::mutate(module_score = meta_df[[module_name]])
    p <- ggplot2::ggplot(embedding_df, ggplot2::aes(x = dim1, y = dim2, color = module_score)) +
         ggplot2::geom_point(size = pt.size,
                             alpha = pt.alpha,
                             stroke = 0) +
         ggplot2::labs(x = embedding_labels[1],
                       y = embedding_labels[2],
                       color = module_label) +
         theme_bayesVG(umap = TRUE)
    if (!is.null(color.palette)) {
      p <- p + ggplot2::scale_color_gradientn(colours = color.palette)
    }
  } else if (plot.type == "violin") {
    p <- ggplot2::ggplot(meta_df, ggplot2::aes(x = !!violin_group_sym, y = !!module_name_sym, color = !!violin_group_sym, fill = !!violin_group_sym))
    if (packageVersion("ggplot2") < "4.0.0") {
      p <- p +
           ggplot2::geom_violin(draw_quantiles = 0.5,
                                scale = "width",
                                alpha = 0.5,
                                linewidth = 0.75)
    } else {
      p <- p +
           ggplot2::geom_violin(quantile.linetype = 1,
                                quantile.linewidth = 0.75,
                                quantiles = 0.5,
                                scale = "width",
                                alpha = 0.5,
                                linewidth = 0.75)
    }
    p <- p +
         ggplot2::labs(x = violin.group,
                       y = module_label,
                       color = violin.group,
                       fill = violin.group) +
         theme_bayesVG()
    if (!is.null(color.palette)) {
      p <- p +
           ggplot2::scale_color_manual(values = color.palette) +
           ggplot2::scale_fill_manual(values = color.palette)
    }
  }
  return(p)
}
