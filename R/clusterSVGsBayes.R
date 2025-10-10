#' Cluster SVGs into modules.
#'
#' @name clusterSVGsBayes
#' @author Jack R. Leary
#' @description This downstream analysis function clusters identified SVGs into modules via an approximate Bayesian soft-clustering approach.
#' @param sp.obj An object of class \code{Seurat} or \code{SpatialExperiment} upon which \code{\link{findSpatiallyVariableFeaturesBayes}} and \code{\link{classifySVGs}} have been run. Defaults to NULL.
#' @param svgs A character vector containing the identified SVGs. Defaults to NULL.
#' @param n.clust An integer specifying the number of clusters to fit to the data. Defaults to 5.
#' @param n.PCs An integer specifying the number of principal components (PCs) to reduce the data to prior to clustering. Defaults to 30.
#' @param algorithm A string specifying the variational inference (VI) approximation algorithm to be used. Must be one of "meanfield", "fullrank", or "pathfinder". Defaults to "fullrank".
#' @param n.iter An integer specifying the maximum number of iterations. Defaults to 30000.
#' @param n.draws An integer specifying the number of draws to be generated from the variational posterior. Defaults to 1000.
#' @param elbo.samples An integer specifying the number of samples to be used to estimate the ELBO at every 100th iteration. Higher values will provide a more accurate estimate at the cost of computational complexity. Defaults to 150 when \code{algorithm} is one of "meanfield" or "fullrank", and 50 when \code{algorithm} is "pathfinder".
#' @param opencl.params A two-element double vector specifying the platform and device IDs of the OpenCL GPU device. Most users should specify \code{c(0, 0)}. See \code{\link[brms]{opencl}} for more details. Defaults to NULL.
#' @param n.cores An integer specifying the number of threads used in compiling and fitting the model and estimating the soft cluster assignment probabilities. Defaults to 2.
#' @param random.seed A double specifying the random seed to be used when fitting and sampling from the model. Defaults to 312.
#' @param verbose A Boolean specifying whether or not verbose model output should be printed to the console. Defaults to TRUE.
#' @details
#' \itemize{
#' \item The soft clustering algorithm is a Gaussian mixture model (GMM), thus each cluster is modeled by a multivariate normal distribution with diagonal covariance. The diagonal covariance assumption is appropriate because each PC is orthogonal to the others.
#' \item The mixing proportions for each cluster are specified with a Dirichlet prior, while the mean follows a Gaussian distribution and the standard deviation a HalfGaussian.
#' \item Due to the architecture of the model, it is necessary for the user to supply a number of clusters via the \code{n.clust} argument. It's difficult to know the correct value to provide beforehand, but luckily the clustering model is quick to run and so multiple values of \code{n.clust} can be fitted and visualized in order to find the "best" value.
#' \item After modeling, the per-cluster assignment probabilities are computed for each gene. The cluster with the highest probability is then defined as the hard cluster for each gene.
#' \item The clustering model also estimates the log-likelihood per-gene, per-draw. This is then used to estimate the overall log-likelihood of the model as well as the corresponding Bayesian information criterion (BIC). The BIC can be used to compare multiple runs of the model, and thus choose the best value of the number of clusters \code{n.clust}.
#' \item When using the fullrank algorithm, it's generally necessary to increase the number of iterations using the \code{n.iter} argument. While the meanfield algorithm generally converges within 1000 iterations, the fullrank algorithm might need e.g., 30,000 iterations. Luckily, the Stan code is very fast so even with 30,000 iterations the clustering should be relatively quick. 
#' }
#' @import magrittr
#' @importFrom cli cli_abort cli_alert_warning
#' @importFrom parallelly availableCores
#' @importFrom SingleCellExperiment logcounts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom irlba prcomp_irlba
#' @importFrom withr with_output_sink
#' @importFrom posterior as_draws_df
#' @importFrom dplyr select mutate rowwise c_across ungroup
#' @importFrom tidyselect starts_with
#' @importFrom matrixStats rowMaxs
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel makeCluster stopCluster
#' @importFrom cli cli_alert_success cli_alert_info
#' @return A list containing the gene-level PCA embedding, a \code{data.frame} of the soft cluster assignments, the fitted model from \code{\link[cmdstanr]{cmdstan_model}}, and the estimated log-likelihood and BIC of the model.
#' @export
#' @examples
#' data(seu_brain)
#' seu_brain <- Seurat::NormalizeData(seu_brain, verbose = FALSE) %>% 
#'              Seurat::FindVariableFeatures(nfeatures = 1000L, verbose = FALSE)
#' seu_brain <- findSpatiallyVariableFeaturesBayes(seu_brain,
#'                                                 naive.hvgs = Seurat::VariableFeatures(seu_brain),
#'                                                 kernel = "matern",
#'                                                 kernel.smoothness = 1.5,
#'                                                 algorithm = "meanfield",
#'                                                 n.cores = 1L,
#'                                                 save.model = TRUE) %>%
#'              classifySVGs(n.SVG = 200L)
#' svg_clusters <- clusterSVGsBayes(seu_brain,
#'                                  svgs = Seurat::VariableFeatures(seu_brain),
#'                                  n.clust = 2L,
#'                                  n.cores = 1L)

clusterSVGsBayes <- function(sp.obj = NULL,
                             svgs = NULL,
                             n.clust = 5L,
                             n.PCs = 30L,
                             algorithm = "fullrank",
                             n.iter = 30000L,
                             n.draws = 1000L,
                             elbo.samples = NULL,
                             opencl.params = NULL,
                             n.cores = 2L,
                             random.seed = 312,
                             verbose = TRUE) {
  # check inputs
  if (is.null(sp.obj) || is.null(svgs)) { cli::cli_abort("All arguments to clusterSVGsBayes() must be supplied.") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (n.clust <= 1L) { cli::cli_abort("Please provide a valid number of clusters.") }
  algorithm <- tolower(algorithm)
  if (!algorithm %in% c("meanfield", "fullrank", "pathfinder")) { cli::cli_abort("Please provide a valid variational inference approximation algorithm.") }
  if (algorithm == "meanfield") { cli::cli_alert_warning("The meanfield VI algorithm often does not perform well on clustering problems do to multi-modality and high correlation of the posterior. Use caution and consider utilizing the default fullrank VI algorithm instead.") }
  if (is.null(elbo.samples)) {
    if (algorithm == "pathfinder") {
      elbo.samples <- 100L
    } else {
      elbo.samples <- 300L
    }
  }
  if (!is.null(opencl.params) && (!is.double(opencl.params) || !length(opencl.params) == 2)) { cli::cli_abort("Argument {.field opencl.params} must be a double vector of length 2 if non-NULL.") }
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
  if (n.cores > unname(parallelly::availableCores())) { cli::cli_abort("The number of requested cores is greater than the number of available cores.") }
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
  # run PCA on the normalized, scaled data
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
                               init = 0.01,
                               num_threads = n.cores,
                               draws = n.draws,
                               opencl_ids = opencl_IDs,
                               num_elbo_draws = elbo.samples,
                               max_lbfgs_iters = 500L,
                               history_size = 100L)
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
                                 init = 0.01,
                                 num_threads = n.cores,
                                 draws = n.draws,
                                 opencl_ids = opencl_IDs,
                                 num_elbo_draws = elbo.samples,
                                 max_lbfgs_iters = 500L,
                                 history_size = 100L)
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
  n_draws <- nrow(theta_draws)
  K <- ncol(theta_draws)
  D <- ncol(mu_draws) / K
  N <- nrow(svg_mtx_pca$x)
  mu_arr <- array(mu_draws, dim = c(n_draws, K, D))
  sigma_arr <- array(sigma_draws, dim = c(n_draws, K, D))
  const_term <- -0.5 * D * log(2 * pi)
  # define function used to compute per-gene, per-draw responsibility values via the log-sum-exp trick
  computeResponsibility <- function(x_i = NULL, 
                                    theta_arr = NULL, 
                                    mu_arr = NULL, 
                                    sigma_arr = NULL) {
    log_probs <- matrix(NA_real_, nrow(theta_arr), ncol(theta_arr))
    for (k in seq_len(ncol(theta_arr))) {
      mu_k <- mu_arr[, k, ]
      sigma_k <- sigma_arr[, k, ]
      diff <- sweep(mu_k, 2, x_i, FUN = "-")
      scaled_sq <- rowSums((diff / sigma_k)^2)
      log_det <- 2 * rowSums(log(sigma_k))
      log_probs[, k] <- log(theta_arr[, k]) - 0.5 * (log_det + scaled_sq)
    }
    log_probs_centered <- log_probs - matrixStats::rowMaxs(log_probs)
    probs <- exp(log_probs_centered)
    probs <- probs / rowSums(probs)
    probs <- colMeans(probs)
    return(probs)
  }
  # compute soft assignment probability for each cluster in parallel
  if (verbose) {
    cli::cli_alert_info("Estimating soft cluster assignment probabilities ...")
  }
  if (verbose) {
    withr::with_output_sink(tempfile(), {
      pb <- utils::txtProgressBar(0, N, style = 3)
    })
    progress_fun <- function(n) utils::setTxtProgressBar(pb, n)
    snow_opts <- list(progress = progress_fun)
  } else {
    snow_opts <- list()
  }
  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl)
  } else {
    cl <- foreach::registerDoSEQ()
  }
  soft_assignments <- foreach(i = seq_len(N),
                              .combine = rbind,
                              .multicombine = ifelse(N > 1, TRUE, FALSE),
                              .maxcombine = ifelse(N > 1, N, 2),
                              .inorder = TRUE,
                              .verbose = FALSE,
                              .options.snow = snow_opts) %dopar% {
    x_i <- svg_mtx_pca$x[i, ]
    computeResponsibilityFast(x_i = x_i, 
                              theta_arr = theta_draws, 
                              mu_arr = mu_arr, 
                              sigma_arr = sigma_arr)
  }
  if (verbose && n.cores > 1L) {
    cat("\n")
  }
  if (n.cores > 1L) {
    parallel::stopCluster(cl)
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
  # estimate log-likelihood from fitted model and summarize
  log_lik_draws <- suppressWarnings(as.data.frame(posterior::as_draws_df(fit_vi$draws(variables = "log_lik")))) %>%
                   dplyr::select(tidyselect::starts_with("log_lik"))
  log_likelihood <- mean(rowSums(log_lik_draws))
  # compute BIC from estimated log-likelihood
  p <- (K - 1) + 2 * K * D
  bic_est <- -2 * log_likelihood + p * log(N)
  if (verbose) {
    cli::cli_alert_success("Posterior summarization complete.")
  }
  # format results
  res <- list(pca_embedding = svg_mtx_pca$x,
              cluster_df = cluster_df,
              model_fit = fit_vi,
              log_likelihood = log_likelihood,
              BIC = bic_est)
  # finish time tracking
  time_diff <- Sys.time() - time_start
  time_units <- ifelse(attributes(time_diff)$units == "secs",
                       "seconds",
                       ifelse(attributes(time_diff)$units == "mins",
                              "minutes",
                              "hours"))
  if (verbose) {
    time_message <- paste0("{.pkg bayesVG} clustering of ",
                           length(svgs),
                           " SVGs completed in ",
                           as.numeric(round(time_diff, 3)),
                           " ",
                           time_units,
                           ".")
    cli::cli_alert_success(time_message)
  }
  return(res)
}
