#' Estimate module score for spatial gene sets. 
#' 
#' @name scoreSpatialModules 
#' @author Jack R. Leary 
#' @description Given a \code{Seurat} or \code{SpatialExperiment} object and the results from \code{\link{clusterSVGsBayes}}, this function uses \code{UCell} to estimate gene set scores for each cluster of SVGs. 
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} upon which \code{\link{findSpatiallyVariableFeaturesBayes}} and \code{\link{classifySVGs}} have been run. Defaults to NULL.
#' @param svg.clusters The results from \code{\link{clusterSVGsBayes}}. Defaults to NULL. 
#' @param n.cores An integer specifying the number of cores used when running \code{UCell}. Defaults to 2. 
#' @importFrom cli cli_abort
#' @importFrom parallelly availableCores
#' @importFrom purrr map
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SpatialExperiment} with per-cluster module scores added to the appropriate metadata slot. Module scores are named as e.g., \code{"svg_cluster_1_UCell"}.
#' @export 
#' @examples
#' data(seu_brain)
#' seu_brain <- Seurat::SCTransform(seu_brain,
#'                                  assay = "Spatial",
#'                                  variable.features.n = 3000L,
#'                                  vst.flavor = "v2",
#'                                  return.only.var.genes = FALSE,
#'                                  seed.use = 312,
#'                                  verbose = FALSE)
#' seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain,
#'                                                 naive.hvgs = Seurat::VariableFeatures(seu_brain),
#'                                                 kernel = "matern",
#'                                                 kernel.smoothness = 1.5,
#'                                                 algorithm = "meanfield",
#'                                                 n.cores = 1L,
#'                                                 save.model = TRUE) %>%
#'              classifySVGs(n.SVG = 1000L)
#' svg_clusters <- clusterSVGsBayes(seu_brain,
#'                                  svgs = Seurat::VariableFeatures(seu_brain),
#'                                  n.clust = 3L,
#'                                  n.cores = 1L)
#' seu_brain <- scoreSpatialModules(seu_brain, 
#'                                  svg.clusters = svg_clusters, 
#'                                  n.cores = 1L)

scoreSpatialModules <- function(sp.obj = NULL, 
                                svg.clusters = NULL, 
                                n.cores = 2L) {
  # check inputs 
  if (is.null(sp.obj) || is.null(svg.clusters)) { cli::cli_abort("Please provide all inputs to scoreSpatialModules().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (n.cores > unname(parallelly::availableCores())) { cli::cli_abort("The number of requested cores is greater than the number of available cores.") }
  # set up list of SVG clusters 
  cluster_list <- split(svg.clusters$cluster_df, svg.clusters$cluster_df$assigned_cluster)
  cluster_list <- purrr::map(cluster_list, \(x) x$gene)
  names(cluster_list) <- paste0("svg_cluster_", seq(length(cluster_list)))
  # run gene set scoring 
  if (inherits(sp.obj, "Seurat")) {
    sp.obj <- UCell::AddModuleScore_UCell(sp.obj, 
                                          features = cluster_list, 
                                          ncores = n.cores)
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    sp.obj <- UCell::ScoreSignatures_UCell(sp.obj, 
                                           features = cluster_list, 
                                           ncores = n.cores)
  }
  return(sp.obj)
}
