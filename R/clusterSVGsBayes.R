#' Cluster SVGs into modules.
#'
#' @name clusterSVGsBayes
#' @author Jack R. Leary
#' @description This downstream analysis function clusters identified SVGs into modules via an approximate Bayesian soft-clustering approach.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} upon which \code{\link{findSpatiallyVariableFeaturesBayes}} and \code{\link{classifySVGs}} have been run. Defaults to NULL.
#' @param svgs A character vector containing the identified SVGs. Defaults to NULL. 
#' @param n.clust An integer specifying the number of clusters to fit to the data. Defaults to 5. 
#' @param n.PCs An integer specifying the number of principal components (PCs) to reduce the data to prior to clustering. Defaults to 30. 
#' @param algorithm A string specifying the variational inference (VI) approximation algorithm to be used. Must be one of "meanfield", "fullrank", or "pathfinder". Defaults to "meanfield".
#' @param n.iter An integer specifying the maximum number of iterations. Defaults to 3000.
#' @param n.draws An integer specifying the number of draws to be generated from the variational posterior. Defaults to 1000.
#' @param elbo.samples An integer specifying the number of samples to be used to estimate the ELBO at every 100th iteration. Higher values will provide a more accurate estimate at the cost of computational complexity. Defaults to 150 when \code{algorithm} is one of "meanfield" or "fullrank", 50 when \code{algorithm} is "pathfinder".  
#' @param opencl.params A two-element double vector specifying the platform and device IDs of the OpenCL GPU device. Most users should specify \code{c(0, 0)}. See \code{\link[brms]{opencl}} for more details. Defaults to NULL.
#' @param n.cores An integer specifying the number of threads used in compiling and fitting the model. Defaults to 2.
#' @param random.seed A double specifying the random seed to be used when fitting and sampling from the model. Defaults to 312.
#' @param verbose A Boolean specifying whether or not verbose model output should be printed to the console. Defaults to TRUE.
#' @details 
#' \itemize{
#' \item The soft clustering algorithm is a Gaussian mixture model (GMM), thus each cluster is modeled by a multivariate normal distribution with diagonal covariance. The diagonal covariance assumption is appropriate because each PC is orthogonal to the others.
#' \item The mixing proportions for each cluster are specified with a Dirichlet prior, while the mean follows a Gaussian distribution and the standard deviation a HalfGaussian. 
#' \item After modeling, the per-cluster assignment probabilities are computed for each gene. The cluster with the highest probability is then defined as the hard cluster for each gene.
#' }
#' @import magrittr 
#' @importFrom parallelly availableCores
#' @importFrom SingleCellExperiment logcounts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom irlba prcomp_irlba
#' @importFrom withr with_output_sink
#' @importFrom posterior as_draws_df
#' @importFrom dplyr select mutate rowwise c_across ungroup
#' @importFrom tidyselect starts_with
#' @importFrom mvtnorm dmvnorm
#' @importFrom cli cli_alert_success
#' @return A \code{data.frame} containing the per-SVG soft cluster assignments.
#' @export
#' @examples
#' data("seu_brain", package = "bayesVG")
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
#'                                  n.clust = 3L)

clusterSVGsBayes <- function(sp.obj = NULL, 
                             svgs = NULL, 
                             n.clust = 5L, 
                             n.PCs = 30L, 
                             algorithm = "meanfield",
                             n.iter = 3000, 
                             n.draws = 1000L, 
                             elbo.samples = NULL, 
                             opencl.params = NULL, 
                             n.cores = 2L, 
                             random.seed = 312, 
                             verbose = TRUE) {
  # check inputs 
  if (is.null(sp.obj) || is.null(svgs)) { stop("All arguments to clusterSVGsBayes() must be supplied.") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { stop("Please provide an object of class Seurat or SpatialExperiment.") }
  algorithm <- tolower(algorithm)
  if (!algorithm %in% c("meanfield", "fullrank", "pathfinder")) { stop("Please provide a valid variational inference approximation algorithm.") }
  if (is.null(elbo.samples)) {
    if (algorithm == "pathfinder") {
      elbo.samples <- 50L
    } else {
      elbo.samples <- 150L
    }
  }
  if (!is.null(opencl.params) && (!is.double(opencl.params) || !length(opencl.params) == 2)) { stop("Argument opencl.params must be a double vector of length 2 if non-NULL.") }
  if (is.null(opencl.params)) {
    opencl_IDs <- NULL
    if (algorithm == "pathfinder") {
      cpp_options <- list(stan_threads = TRUE)
    } else {
      cpp_options <- list(stan_threads = FALSE)
    }
  } else {
    opencl_IDs <- opencl.params
    if (algorithm == "pathfinder") {
      cpp_options <- list(stan_opencl = TRUE, stan_threads = TRUE)
    } else {
      cpp_options <- list(stan_opencl = TRUE, stan_threads = FALSE)
    }
  }
  if (n.cores > unname(parallelly::availableCores())) { stop("The number of requested cores is greater than the number of available cores.") }
  # start time tracking 
  time_start <- Sys.time()
  # extract matrix of normalized gene expression & scale it 
  if (inherits(sp.obj, "Seurat")) {
    expr_mtx <- Seurat::GetAssayData(sp.obj,
                                     assay = Seurat::DefaultAssay(sp.obj),
                                     layer = "data")
  } else {
    expr_mtx <- SingleCellExperiment::logcounts(sp.obj)
  }
  expr_mtx <- as.matrix(expr_mtx[svgs, ])
  expr_mtx <- t(scale(t(expr_mtx)))
  attributes(expr_mtx)[3:4] <- NULL
  # run PCA on the normalized data 
  svg_mtx_pca <- irlba::prcomp_irlba(expr_mtx, 
                                     n = n.PCs, 
                                     center = FALSE, 
                                     scale. = FALSE)
  # prepare data to be passed to cmdstan
  data_list <- list(N = nrow(svg_mtx_pca$x), 
                    D = ncol(svg_mtx_pca$x), 
                    K = n.clust, 
                    X = svg_mtx_pca$x)
  stan_file <- system.file("clusterSVGs.stan", package = "bayesVG")
  # compile model
  mod <- cmdstan_model(stan_file, compile = FALSE)
  mod$compile(cpp_options = cpp_options,
              stanc_options = list("O1"),
              force_recompile = TRUE,
              threads = n.cores)
  # fit model with desired algorithm
  if (verbose) {
    if (algorithm %in% c("meanfield", "fullrank")) {
      fit_vi <- mod$variational(data_list,
                                seed = random.seed,
                                init = 0,
                                algorithm = algorithm,
                                iter =  n.iter,
                                draws = n.draws,
                                opencl_ids = opencl_IDs,
                                elbo_samples = elbo.samples)
    } else {
      fit_vi <- mod$pathfinder(data_list,
                               seed = random.seed,
                               init = 0, 
                               num_threads = n.cores,
                               draws = n.draws,
                               opencl_ids = opencl_IDs,
                               num_elbo_draws = elbo.samples, 
                               max_lbfgs_iters = 100L, 
                               history_size = 25L)
    }
  } else {
    withr::with_output_sink(tempfile(), {
      if (algorithm %in% c("meanfield", "fullrank")) {
        fit_vi <- mod$variational(data_list,
                                  seed = random.seed,
                                  init = 0,
                                  algorithm = algorithm,
                                  iter =  n.iter,
                                  draws = n.draws,
                                  opencl_ids = opencl_IDs,
                                  elbo_samples = elbo.samples)
      } else {
        fit_vi <- mod$pathfinder(data_list,
                                 seed = random.seed,
                                 init = 0, 
                                 num_threads = n.cores,
                                 draws = n.draws,
                                 opencl_ids = opencl_IDs,
                                 num_elbo_draws = elbo.samples, 
                                 max_lbfgs_iters = 100L, 
                                 history_size = 25L)
      }
    })
  }
  # optionally run diagnostics
  if (verbose) {
    fit_vi$cmdstan_diagnose()
  }
  # extract posterior draws 
  draws_df <- posterior::as_draws_df(fit_vi$draws())
  theta_draws <- suppressWarnings(as.matrix(dplyr::select(draws_df, tidyselect::starts_with("theta"))))
  mu_draws    <- suppressWarnings(as.matrix(dplyr::select(draws_df, tidyselect::starts_with("mu"))))
  sigma_draws <- suppressWarnings(as.matrix(dplyr::select(draws_df, tidyselect::starts_with("sigma"))))
  # set up constants
  n_iter <- nrow(theta_draws)
  K <- ncol(theta_draws)
  D <- ncol(mu_draws) / K
  N <- nrow(svg_mtx_pca$x)
  # define function used to compute per-gene, per-draw responsibility values
  computeResponsibility <- function(x = NULL, 
                                    theta_s = NULL, 
                                    mu_s = NULL, 
                                    sigma_s = NULL,
                                    K = NULL, 
                                    D = NULL) {
    log_probs <- vector("numeric", length = K)
    for (i in seq(K)) {
      idx <- ((i - 1) * D + 1):(i * D)
      mu_i <- mu_s[idx]
      sigma_i <- sigma_s[idx]
      cov_i <- diag(sigma_i^2)
      log_probs[i] <- log(theta_s[i]) + mvtnorm::dmvnorm(x, mean = mu_i, sigma = cov_i, log = TRUE)
    }
    max_log <- max(log_probs)
    probs <- exp(log_probs - (max_log + log(sum(exp(log_probs - max_log)))))
    return(probs)
  }
  # compute soft assignment probability for each cluster 
  soft_assignments <- matrix(0, nrow = N, ncol = K)
  for (i in seq(N)) {
    x_i <- as.numeric(svg_mtx_pca$x[i, ])
    resp_iter <- matrix(NA, nrow = n_iter, ncol = K)
    for (s in seq(n_iter)) {
      theta_s <- theta_draws[s, ]
      mu_s <- mu_draws[s, ]
      sigma_s <- sigma_draws[s, ]
      resp_iter[s, ] <- computeResponsibility(x_i, 
                                              theta_s = theta_s, 
                                              mu_s = mu_s, 
                                              sigma_s = sigma_s, 
                                              K = K, 
                                              D = D)
    }
    soft_assignments[i, ] <- colMeans(resp_iter)
  }
  # generate a hard assignment for each gene to its most likely cluster
  cluster_df <- as.data.frame(soft_assignments)
  colnames(cluster_df) <- paste0("prob_cluster_", seq(K))
  cluster_df <- dplyr::mutate(cluster_df, 
                              gene = svgs, 
                              .before = 1) %>% 
                dplyr::rowwise() %>%
                dplyr::mutate(assigned_cluster = which.max(dplyr::c_across(tidyselect::starts_with("prob_cluster_")))) %>%
                dplyr::ungroup()
  if (verbose) {
    cli::cli_alert_success("Posterior summarization complete.")
  }
  # finish time tracking 
  time_diff <- Sys.time() - time_start
  time_units <- ifelse(attributes(time_diff)$units == "secs", 
                       "seconds", 
                       ifelse(attributes(time_diff)$units == "mins", 
                              "minutes", 
                              "hours"))
  if (verbose) {
    time_message <- paste0("bayesVG clustering of ", 
                           length(svgs), 
                           " SVGs completed in ", 
                           as.numeric(round(time_diff, 3)), 
                           " ", 
                           time_units, 
                           ".")
    cli::cli_alert_success(time_message)
  }
  return(cluster_df)
}
