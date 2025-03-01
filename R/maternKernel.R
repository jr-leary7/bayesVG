#' Compute the Matern-family kernel.
#'
#' @name maternKernel
#' @author Jack R. Leary
#' @description This function computes the Matern kernel for a given vector of distances and a global length-scale.
#' @param d A numeric vector of distances. Defaults to NULL. 
#' @param length.scale A numeric specifying the length-scale \eqn{\ell} for the approximate GP. Defaults to NULL. 
#' @param nu A numeric specifying the smoothness parameter \eqn{\nu} of the Matern kernel. Defaults to NULL. 
#' @param sigma A numeric specifying the variance of the Matern kernel. Defaults to 1. 
#' @return A vector of numeric values. 
#' @seealso \code{\link{expQuadKernel}}
#' @seealso \code{\link{periodicKernel}}

maternKernel <- function(d = NULL,
                         length.scale = NULL,
                         nu = NULL,
                         sigma = 1) {
  # check inputs 
  if (is.null(d) || is.null(length.scale) || is.null(nu)) { stop("All parameters must be provided.") }
  # compute scaled distances
  scaled_distance <- sqrt(2 * nu) * d / length.scale
  # compute Matern kernel values
  res <- sigma^2 * (2^(1 - nu) / gamma(nu)) * (scaled_distance)^nu * besselK(scaled_distance, nu)
  # replace NA values with sigma^2
  res[is.nan(res)] <- sigma^2
  return(res)
}
