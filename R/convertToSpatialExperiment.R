#' Convert from \code{Seurat} to \code{SpatialExperiment}.
#'
#' @name convertToSpatialExperiment
#' @author Jack R. Leary
#' @description This function converts a \code{Seurat} object to a \code{SpatialExperiment} object while retaining all assays, metadata, and images.
#' @param seu.obj An object of class \code{Seurat}. Defaults to NULL.
#' @param sample.id A string specifying the sample ID corresponding to the image saved in \code{seu.obj}.
#' @param scale.coords A Boolean specifying whether the spatial coordinates matrix should be scaled. Defaults to FALSE. 
#' @importFrom cli cli_abort
#' @importFrom Seurat as.SingleCellExperiment GetTissueCoordinates
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom grDevices as.raster
#' @importFrom SpatialExperiment SpatialImage SpatialExperiment
#' @importFrom SummarizedExperiment assays rowData colData
#' @importFrom SingleCellExperiment reducedDims altExps
#' @importFrom dplyr select
#' @importFrom coop scaler
#' @return An object of class \code{SpatialExperiment}.
#' @seealso \code{\link[Seurat]{as.SingleCellExperiment}}
#' @seealso \code{\link[SpatialExperiment]{SpatialExperiment}}
#' @export
#' @examples
#' data(seu_brain)
#' spe_brain <- convertToSpatialExperiment(seu_brain, sample.id = "anterior1")

convertToSpatialExperiment <- function(seu.obj = NULL,
                                       sample.id = NULL,
                                       scale.coords = FALSE) {
  # check inputs
  if (!inherits(seu.obj, "Seurat")) { cli::cli_abort("Please provide an object of class {.pkg Seurat}.") }
  if (!sample.id %in% names(seu.obj@images)) { cli::cli_abort("Please provide a valid {.field sample.id} value.") }
  # convert to singlecellexperiment object
  sce <- Seurat::as.SingleCellExperiment(seu.obj)
  # extract image data
  img_data <- S4Vectors::DataFrame(sample_id = sample.id,
                                   image_id = sample.id,
                                   data = I(list(SpatialExperiment::SpatialImage(grDevices::as.raster(seu.obj@images[[sample.id]]@image)))),
                                   scaleFactor = seu.obj@images[[sample.id]]@scale.factors$lowres)
  # extract spatial coordinates and optionally scale them
  spatial_coords <- as.matrix(dplyr::select(Seurat::GetTissueCoordinates(seu.obj), 1:2))
  if (scale.coords) {
    spatial_coords <- coop::scaler(spatial_coords)
    attributes(spatial_coords)[3:4] <- NULL
  }
  # make sure to retain correct dimnames
  colnames(spatial_coords) <- c("x", "y")
  rownames(spatial_coords) <- colnames(seu.obj)
  # create spatialexperiment object
  spe <- SpatialExperiment::SpatialExperiment(assays = SummarizedExperiment::assays(sce),
                                              rowData = SummarizedExperiment::rowData(sce),
                                              colData = SummarizedExperiment::colData(sce),
                                              metadata = S4Vectors::metadata(sce),
                                              reducedDims = SingleCellExperiment::reducedDims(sce),
                                              altExps = SingleCellExperiment::altExps(sce),
                                              sample_id = sample.id,
                                              spatialCoords = spatial_coords,
                                              imgData = img_data)
  spe$in_tissue <- 1
  spe$sample_id <- sample.id
  return(spe)
}
