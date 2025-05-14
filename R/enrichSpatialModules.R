#' Run enrichment analysis on SVG modules. 
#' 
#' @name enrichSpatialModules
#' @author Jack R. Leary
#' @description This function utilizes \code{\link[gprofiler2]{gost}} to perform enrichment analysis on each SVG module. 
#' @param svg.clusters The results from \code{\link{clusterSVGsBayes}}. Defaults to NULL.
#' @param species A string specifying the species from which the cells originate. Defaults to "hsapiens". 
#' @importFrom cli cli_abort
#' @importFrom purrr map
#' @importFrom rlang sym
#' @importFrom dplyr arrange desc 
#' @return A \code{data.frame} with enrichment analysis results for each spatial module. 
#' @export 
#' @examples
#' data(seu_brain)
#' seu_brain <- Seurat::NormalizeData(seu_brain, verbose = FALSE) %>% 
#'              Seurat::FindVariableFeatures(nfeatures = 3000L, verbose = FALSE)
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
#' enrich_res <- enrichSpatialModules(svg_clusters, species = "mmusculus")

enrichSpatialModules <- function(svg.clusters = NULL, species = "hsapiens") {
  # check inputs 
  if (is.null(svg.clusters)) { cli::cli_abort("Please provide all inputs to enrichSpatialModules().") }
  # run enrichment analysis 
  module_list <- split(svg.clusters$cluster_df, svg.clusters$cluster_df$assigned_cluster)
  module_list <- purrr::map(module_list, \(x) {
    cluster_ID <- unique(x$assigned_cluster)
    correct_prob_column <- paste0("prob_cluster_", cluster_ID)
    correct_prob_column_sym <- rlang::sym(correct_prob_column)
    x <- dplyr::arrange(x, dplyr::desc(!!correct_prob_column_sym))
    return(x$gene)
  })
  enrich_res <- gprofiler2::gost(module_list, 
                                 organism = species, 
                                 ordered_query = TRUE, 
                                 multi_query = FALSE, 
                                 significant = FALSE)$result
  return(enrich_res)
}
