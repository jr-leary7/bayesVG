#' Compute the periodic kernel.
#'
#' @name periodicKernel
#' @author Jack R. Leary
#' @description This function computes the exponentiated quadratic kernel for a given vector of distances and a global length-scale.
#' @param d A numeric vector of distances. Defaults to NULL. 
#' @param length.scale A numeric specifying the length-scale \eqn{\ell} for the approximate GP. Defaults to NULL. 
#' @param period An integer specifying the period of the kernel. Defaults to 100.  
#' @return A vector of numeric values. 
#' @seealso \code{\link{maternKernel}}
#' @seealso \code{\link{expQuadKernel}}

periodicKernel <- function(d = NULL, 
                           length.scale = NULL,
                           period = 100L) {
  # check inputs 
  if (is.null(d) || is.null(length.scale)) { stop("All parameters must be provided.") }
  # compute periodic kernel values
  d_sqrt <- sqrt(d)
  res <- exp(-2 * sin(pi * d_sqrt / period)^2 / (length.scale^2))
  return(res)
}
