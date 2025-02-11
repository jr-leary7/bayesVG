#' Extract the fitted model.
#'
#' @name extractModel
#' @author Jack R. Leary
#' @description This function extracts the fitted model from \code{brms} or \code{cmdstanr} (depending on whether you're identifying HVGs or SVGs) from the user-provided object to which the model has been written.
#' @param obj The \code{Seurat}, \code{SingleCellExperiment}, or \code{SpatialExperiment} whose unstructured metadata contains the fitted model. Defaults to NULL.
#' @importFrom Seurat DefaultAssay
#' @details
#' \itemize{
#' \item In order for the fitted model to have been saved, the modeling function must have been called with the argument \code{save.model = TRUE}. Otherwise, the fitted model is discarded from memory after the modeling function has finished running.
#' }
#' @seealso \code{findVariableFeaturesBayes}
#' @seealso \code{findSpatiallyVariableFeaturesBayes}
#' @export

extractModel <- function(obj = NULL) {
  # check inputs
  if (is.null(obj)) { stop("Argument obj must be non-NULL.") }
  # extract fitted model from object's unstructured metadata
  if (inherits(obj, "Seurat")) {
    model_fit <- obj@assays[[Seurat::DefaultAssay(obj)]]@misc$model_fit
  } else if (inherits(obj, "SingleCellExperiment") || inherits(obj, "SpatialExperiment")) {
    model_fit <- obj@metadata$model_fit
  }
  return(model_fit)
}
