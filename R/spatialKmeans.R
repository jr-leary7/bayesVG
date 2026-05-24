#' Run \emph{k}-means clustering on spatial coordinates with special conditions.
#'
#' @name spatialKmeans
#' @author Jack R. Leary
#' @description This helper function runs \emph{k}-means clustering on the spatial coordinates matrix, with added functionalities that make convergence more likely.
#' @param coord_mtx A matrix of spatial coordinates. Defaults to NULL. 
#' @param iter.max An integer specifying the maximum number of algorithm iterations. Defaults to 100. 
#' @param n.start An integer passed to the `nstart` argument of \code{\link[stats]{kmeans}}. Defaults to 10. 
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom stats kmeans 
#' @return An object of class \code{kmeans}.
#' @seealso \code{\link[stats]{kmeans}}

spatialKmeans <- function(coord.mtx = NULL, 
                          n.centers = 20L, 
                          iter.max = 100L, 
                          n.start = 10L) {
  # check inputs
  if (is.null(coord.mtx) || !inherits(coord.mtx, "matrix")) { cli::cli_abort("You must provide a matrix of spatial coordinates to spatialKmeans().") }
  # run hartigan-wong k-means with fallback to (slower) lloyd algorithm on qt warning 
  evaluation_env <- new.env(parent = emptyenv())
  evaluation_env$qt_warning_thrown <- FALSE
  kmeans_res <- withCallingHandlers(
    stats::kmeans(coord.mtx, 
                  centers = n.centers, 
                  iter.max = iter.max, 
                  nstart = n.start, 
                  algorithm = "Hartigan-Wong"
    ),
    warning = \(w) {
      if (grepl("Quick-TRANSfer", w$message)) {
        evaluation_env$qt_warning_thrown <- TRUE
        if (!is.null(findRestart("muffleWarning"))) {
          invokeRestart("muffleWarning")
        }
      }
    }
  )
  if (evaluation_env$qt_warning_thrown) {
    cli::cli_alert_info("The Hartigan-Wong k-means algorithm hit its QT limit. Falling back to Lloyd's k-means algorithm...")
    kmeans_res <- stats::kmeans(coord.mtx, 
                                centers = n.centers, 
                                iter.max = iter.max, 
                                nstart = n.start, 
                                algorithm = "Lloyd")
  }
  return(kmeans_res)
}
