#' Convert from \code{Seurat} to \code{SpatialExperiment}.
#'
#' @name convertToSpatialExperiment
#' @author Jack R. Leary
#' @description This function converts a \code{Seurat} object to a \code{SpatialExperiment} object while retaining all assays, metadata, and images.
#' @param seu.obj An object of class \code{Seurat}. Defaults to NULL.
#' @param sample.id A string specifying the sample ID corresponding to the image saved in \code{seu.obj}.
#' @param scale.coords A Boolean specifying whether the spatial coordinates matrix should be scaled. Defaults to FALSE.
#' @importFrom Seurat as.SingleCellExperiment GetTissueCoordinates
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SpatialExperiment SpatialImage SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom SingleCellExperiment reducedDims altExps
#' @importFrom dplyr select
#' @return An object of class \code{SpatialExperiment}.
#' @seealso \code{\link[Seurat]{as.SingleCellExperiment}}
#' @seealso \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @export

convertToSpatialExperiment <- function(seu.obj = NULL,
                                       sample.id = NULL,
                                       scale.coords = FALSE) {
  # check inputs
  if (!inherits(seu.obj, "Seurat")) { stop("Please provide an object of class Seurat.") }
  if (!sample.id %in% names(seu.obj@images)) { stop("Please provide a valid sample.id value.") }
  # convert to singlecellexperiment object
  sce <- Seurat::as.SingleCellExperiment(seu.obj)
  # extract image data
  img_data <- S4Vectors::DataFrame(sample_id = sample.id,
                                   image_id = sample.id,
                                   data = I(list(SpatialExperiment::SpatialImage(as.raster(seu.obj@images[[sample.id]]@image)))),
                                   scaleFactor = seu.obj@images[[sample.id]]@scale.factors$lowres)
  # extract spatial coordinates and optionally scale them
  spatial_coords <- as.matrix(dplyr::select(Seurat::GetTissueCoordinates(seu.obj), -cell))
  if (scale.coords) {
    spatial_coords <- scale(spatial_coords)
  }
  # create spatialexperiment object
  spe <- SpatialExperiment::SpatialExperiment(assays = SummarizedExperiment::assays(sce_brain),
                                              rowData = SummarizedExperiment::rowData(sce_brain),
                                              colData = SummarizedExperiment::colData(sce_brain),
                                              metadata = S4vectors::metadata(sce_brain),
                                              reducedDims = SingleCellExperiment::reducedDims(sce_brain),
                                              altExps = SingleCellExperiment::altExps(sce_brain),
                                              sample_id = sample.id,
                                              spatialCoords = spatial_coords,
                                              imgData = img_data)
  spe$in_tissue <- 1
  spe$sample_id <- sample.id
  return(spe)
}
