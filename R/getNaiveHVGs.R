#' Compute a set of naive HVGs. 
#' 
#' @name getNaiveHVGs
#' @author Jack R. Leary
#' @description This function, given an object of class \code{Seurat} or \code{SpatialExperiment}, computes a set of naive HVGs to be used as candidates for SVG modeling.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} containing spatial data. The function assumes that the raw counts have already been depth- and log1p-normalized. Defaults to NULL.
#' @param n.hvg An integer specifying the number of naive HVGs to be returned. Defaults to 3000. 
#' @importFrom cli cli_abort
#' @importFrom Seurat FindVariableFeatures VariableFeatures
#' @importFrom scran modelGeneVar getTopHVGs
#' @return A character vector of length \code{n.hvg} specifying the set of naive HVGs. 
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @seealso \code{\link[scran]{modelGeneVar}}
#' @export 
#' @examples 
#' data(seu_brain)
#' seu_brain <- Seurat::NormalizeData(seu_brain, verbose = FALSE)
#' naive_hvgs <- getNaiveHVGs(seu_brain)

getNaiveHVGs <- function(sp.obj = NULL, n.hvg = 3000L) {
  # check inputs 
  if (is.null(sp.obj)) { cli::cli_abort("Please provide a spatial data object to getNaiveHVGs().") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  # compute set of naive HVGs
  if (inherits(sp.obj, "Seurat")) {
    sp.obj <- Seurat::FindVariableFeatures(sp.obj, 
                                           nfeatures = n.hvg, 
                                           verbose = FALSE)
    naive_hvgs <- Seurat::VariableFeatures(sp.obj)
  } else if (inherits(sp.obj, "SpatialExperiment")) {
    naive_hvgs <- scran::getTopHVGs(scran::modelGeneVar(sp.obj), n = n.hvg)
  }
  return(naive_hvgs)
}
