#' Compute the exponentiated quadratic kernel.
#'
#' @name expQuadKernel
#' @author Jack R. Leary
#' @description This function computes the exponentiated quadratic kernel for a given vector of distances and a global length-scale.
#' @param d A numeric vector of distances. Defaults to NULL. 
#' @param length.scale A numeric specifying the length-scale \eqn{\ell} for the approximate GP. Defaults to NULL. 
#' @return A vector of numeric values. 
#' @seealso \code{\link{maternKernel}}

expQuadKernel <- function(d = NULL, length.scale = NULL) {
  # check inputs 
  if (is.null(d) || is.null(length.scale)) { stop("All parameters must be provided.") }
  # compute exponentiated quadratic kernel values
  res <- exp(-d / (2 * length.scale^2))
  return(res)
}
