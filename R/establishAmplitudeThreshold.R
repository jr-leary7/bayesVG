#' Identify spatially variable genes in a Bayesian manner.
#'
#' @name establishAmplitudeThreshold
#' @author Jack R. Leary
#' @description This function fits the same Gaussian process model as \code{\link{findSpatiallyVariableFeaturesBayes}}, but does so for a subset of genes on shuffled spatial coordinates.
#' @param sp.obj An object of class \code{Seurat} containing spatial data. Defaults to NULL.
#' @param sampled.genes A vector containing a subsample of the rownames of \code{sp.obj}. Defaults to NULL.
#' @param kernel A string specifying the covariance kernel to be used when fitting the GP. Must be one of "exp_quad", "matern", or "periodic". Defaults to "exp_quad".
#' @param kernel.smoothness A double specifying the smoothness parameter \eqn{\nu} used when computing the Matérn kernel. Must be one of 0.5, 1.5, or 2.5. Using 0.5 corresponds to the exponential kernel. Defaults to 1.5.
#' @param kernel.period An integer specifying the period parameter \eqn{p} used when computing the periodic kernel. Defaults to 100.
#' @param n.basis.fns An integer specifying the number of basis functions to be used when approximating the GP as a Hilbert space. Defaults to 20.
#' @param algorithm A string specifying the variational inference (VI) approximation algorithm to be used. Must be one of "meanfield", "fullrank", or "pathfinder". Defaults to "meanfield".
#' @param mle.init A Boolean specifying whether the the VI algorithm should be initialized using the MLE for each parameter. In general, this isn't strictly necessary but can help if the VI algorithm struggles to converge when provided with the default initialization (zero). Cannot be used when the Pathfinder algorithm is specified. Defaults to TRUE.
#' @param gene.depth.adjust A Boolean specifying whether the model should include a fixed effect term for total gene expression. Defaults to TRUE.
#' @param n.cores An integer specifying the number of threads used in compiling and fitting the model. Defaults to 2.
#' @param random.seed A double specifying the random seed to be used when fitting and sampling from the model. Defaults to 312.
#' @param verbose A Boolean specifying whether or not verbose model output should be printed to the console. Defaults to TRUE.
#' @details
#' \itemize{
#' \item This function is meant to be run in tandem with \code{\link{findSpatiallyVariableFeaturesBayes}} in order to determine an accurate threshold for the gene-specific amplitude parameter \eqn{\tau_g}.
#' }
#' @import magrittr
#' @import cmdstanr
#' @importFrom cli cli_abort cli_alert_warning cli_alert_success
#' @importFrom parallelly availableCores
#' @importFrom Seurat GetAssayData DefaultAssay GetTissueCoordinates
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom BiocGenerics counts
#' @importFrom Matrix rowSums
#' @importFrom dplyr relocate mutate rename rename_with with_groups select inner_join desc filter distinct arrange
#' @importFrom tidyr pivot_longer
#' @importFrom stats kmeans dist median
#' @importFrom withr with_output_sink
#' @return A \code{data.frame} containing amplitude-related statistics for each gene.
#' @seealso \code{\link{findSpatiallyVariableFeaturesBayes}}
#' @export
#' @examples
#' data(seu_brain)
#' seu_brain <- Seurat::SCTransform(seu_brain,
#'                                  assay = "Spatial",
#'                                  variable.features.n = 3000L,
#'                                  vst.flavor = "v2",
#'                                  return.only.var.genes = FALSE,
#'                                  seed.use = 312,
#'                                  verbose = FALSE)
#' set.seed(312)
#' sampled_genes <- sample(VariableFeatures(seu_brain), size = 50L)
#' amplitude_summary <- establishAmplitudeThreshold(seu_brain, sampled.genes = sampled_genes)

establishAmplitudeThreshold <- function(sp.obj = NULL,
                                        sampled.genes = NULL,
                                        kernel = "exp_quad",
                                        kernel.smoothness = 1.5,
                                        kernel.period = 100L,
                                        n.basis.fns = 20L,
                                        algorithm = "meanfield",
                                        mle.init = TRUE,
                                        gene.depth.adjust = TRUE,
                                        n.cores = 2L,
                                        random.seed = 312,
                                        verbose = TRUE) {
  # check inputs
  if (is.null(sp.obj)) { cli::cli_abort("{.field sp.obj} must be non-NULL.") }
  if (!(inherits(sp.obj, "Seurat") || inherits(sp.obj, "SpatialExperiment"))) { cli::cli_abort("Please provide an object of class {.pkg Seurat} or {.pkg SpatialExperiment}.") }
  if (is.null(sampled.genes)) {
    cli::cli_alert_info("Sampling 50 genes from the rownames of {.field sp.obj} ...")
    sampled.genes <- sample(rownames(sp.obj),
                            size = 50L,
                            replace = FALSE)
  }
  kernel <- tolower(kernel)
  if (!kernel %in% c("exp_quad", "matern", "periodic")) { cli::cli_abort("Please provide a valid covariance kernel.") }
  if (kernel == "matern" && !kernel.smoothness %in% c(0.5, 1.5, 2.5)) { cli::cli_abort("When utilizing the Matérn kernel you must provide a valid smoothness parameter value.") }
  algorithm <- tolower(algorithm)
  if (!algorithm %in% c("meanfield", "fullrank", "pathfinder")) { cli::cli_abort("Please provide a valid variational inference approximation algorithm.") }
  if (mle.init && algorithm == "pathfinder") { cli::cli_alert_warning("Initialization at the MLE is not supported when using the Pathfinder algorithm.") }
  if (n.cores > unname(parallelly::availableCores())) { cli::cli_abort("The number of requested cores is greater than the number of available cores.") }
  # set up ELBO sampling strategy & C++ options
  if (algorithm == "pathfinder") {
    elbo_samples <- 50L
    cpp_options <- list(stan_threads = TRUE)
  } else {
    elbo_samples <- 150L
    cpp_options <- list(stan_threads = FALSE)
  }
  # extract spatial coordinates & scale them
  if (inherits(sp.obj, "Seurat")) {
    spatial_df <- Seurat::GetTissueCoordinates(sp.obj) %>%
                  dplyr::select(1:2)
  } else {
    spatial_df <- SpatialExperiment::spatialCoords(sp.obj) %>%
                  as.data.frame()
  }
  spatial_mtx <- scale(as.matrix(spatial_df))
  # shuffle spatial coordinates
  spatial_mtx <- spatial_mtx[sample(nrow(spatial_mtx)), ]
  # extract matrix of normalized gene expression
  if (inherits(sp.obj, "Seurat")) {
    expr_mtx <- Seurat::GetAssayData(sp.obj,
                                     assay = Seurat::DefaultAssay(sp.obj),
                                     layer = "data")
  } else {
    expr_mtx <- SingleCellExperiment::logcounts(sp.obj)
  }
  # convert expression matrix to long data.frame for modeling & post-process
  expr_df <- as.data.frame(expr_mtx[sampled.genes, ]) %>%
             dplyr::mutate(gene = rownames(.), .before = 1) %>%
             tidyr::pivot_longer(cols = !gene,
                                 names_to = "spot",
                                 values_to = "gene_expression") %>%
             dplyr::relocate(spot, gene) %>%
             dplyr::mutate(gene = factor(gene, levels = unique(gene)),
                           spot = factor(spot, levels = unique(spot)),
                           gene_expression = as.numeric(scale(gene_expression))) %>%
             as.data.frame()
  # estimate global length-scale
  M <- nrow(spatial_mtx)
  kmeans_centers <- stats::kmeans(spatial_mtx, centers = n.basis.fns, iter.max = 100L)$centers
  dists_centers <- as.matrix(stats::dist(kmeans_centers))
  lscale <- stats::median(dists_centers[upper.tri(dists_centers)])
  # estimate matrix of basis functions used to approximate GP with desired kernel
  phi <- matrix(0, nrow = M, ncol = n.basis.fns)
  for (i in seq(n.basis.fns)) {
    dist_vec <- rowSums((spatial_mtx - matrix(kmeans_centers[i, ], nrow = M, ncol = 2, byrow = TRUE))^2)
    if (kernel == "exp_quad") {
      phi[, i] <- expQuadKernel(dist_vec, length.scale = lscale)
    } else if (kernel == "matern") {
      phi[, i] <- maternKernel(dist_vec,
                               length.scale = lscale,
                               nu = kernel.smoothness)
    } else if (kernel == "periodic") {
      phi[, i] <- periodicKernel(dist_vec,
                                 length.scale = lscale,
                                 period = kernel.period)
    }
  }
  # scale basis functions
  phi <- scale(phi)
  attributes(phi)[2:3] <- NULL
  # compute some constants
  N <- nrow(expr_df)
  G <- length(unique(expr_df$gene))
  # prepare data to be passed to cmdstan
  if (gene.depth.adjust) {
    if (inherits(sp.obj, "Seurat")) {
      expr_tmp <- Seurat::GetAssayData(sp.obj,
                                       assay = Seurat::DefaultAssay(sp.obj),
                                       layer = "counts")
    } else {
      expr_tmp <- BiocGenerics::counts(sp.obj)
    }
    expr_tmp <- expr_tmp[sampled.genes, ]
    gene_depths <- unname(log(Matrix::rowSums(expr_tmp)))
    data_list <- list(M = M,
                      N = N,
                      G = G,
                      k = n.basis.fns,
                      spot_id = as.integer(expr_df$spot),
                      gene_id = as.integer(expr_df$gene),
                      phi = phi,
                      gene_depths = gene_depths,
                      y = expr_df$gene_expression)
    stan_file <- system.file("approxGP2.stan", package = "bayesVG")
  } else {
    data_list <- list(M = M,
                      N = N,
                      G = G,
                      k = n.basis.fns,
                      spot_id = as.integer(expr_df$spot),
                      gene_id = as.integer(expr_df$gene),
                      phi = phi,
                      y = expr_df$gene_expression)
    stan_file <- system.file("approxGP.stan", package = "bayesVG")
  }
  # compile model
  mod <- cmdstan_model(stan_file, compile = FALSE)
  mod$compile(cpp_options = cpp_options,
              stanc_options = list("O1"),
              force_recompile = TRUE,
              threads = n.cores)
  # fit model with desired algorithm
  if (verbose) {
    if (algorithm %in% c("meanfield", "fullrank")) {
      if (mle.init) {
        model_init <- mod$optimize(data_list,
                                   seed = random.seed,
                                   init = 0,
                                   jacobian = FALSE,
                                   iter = 1000L,
                                   algorithm = "lbfgs",
                                   history_size = 25L)
      } else {
        model_init <- 0
      }
      fit_vi <- mod$variational(data_list,
                                seed = random.seed,
                                init = model_init,
                                algorithm = algorithm,
                                iter =  3000L,
                                draws = 1000L,
                                elbo_samples = elbo_samples)
    } else {
      fit_vi <- mod$pathfinder(data_list,
                               seed = random.seed,
                               init = 0.1,
                               num_threads = n.cores,
                               draws = 100L,
                               num_paths = n.cores,
                               max_lbfgs_iters = 200L,
                               num_elbo_draws = elbo_samples,
                               history_size = 25L)
    }
  } else {
    withr::with_output_sink(tempfile(), {
      if (algorithm %in% c("meanfield", "fullrank")) {
        if (mle.init) {
          model_init <- mod$optimize(data_list,
                                     seed = random.seed,
                                     init = 0,
                                     jacobian = FALSE,
                                     iter = 1000L,
                                     algorithm = "lbfgs",
                                     history_size = 25L)
        } else {
          model_init <- 0
        }
        fit_vi <- mod$variational(data_list,
                                  seed = random.seed,
                                  init = model_init,
                                  algorithm = algorithm,
                                  iter =  3000L,
                                  draws = 1000L,
                                  elbo_samples = elbo_samples)
      } else {
        fit_vi <- mod$pathfinder(data_list,
                                 seed = random.seed,
                                 init = 0.1,
                                 num_threads = n.cores,
                                 draws = 1000L,
                                 num_paths = n.cores,
                                 max_lbfgs_iters = 200L,
                                 num_elbo_draws = elbo_samples,
                                 history_size = 25L)
      }
    })
  }
  # summarize posterior
  gene_mapping <- data.frame(gene = as.character(expr_df$gene),
                             gene_id = as.character(as.integer(expr_df$gene))) %>%
                  dplyr::distinct()
  amplitude_summary <- fit_vi$summary(variables = "amplitude") %>%
                       dplyr::rename_with(~paste0("amplitude_", .), .cols = -1) %>%
                       dplyr::rename(amplitude_ci_ll = amplitude_q5,
                                     amplitude_ci_ul = amplitude_q95) %>%
                       dplyr::mutate(gene_id = sub("^.*\\[(.*)\\].*$", "\\1", variable), .before = 1) %>%
                       dplyr::inner_join(gene_mapping, by = "gene_id") %>%
                       dplyr::relocate(gene) %>%
                       dplyr::select(-c(variable, gene_id)) %>%
                       dplyr::mutate(amplitude_dispersion = amplitude_sd^2 / amplitude_mean) %>%
                       dplyr::arrange(dplyr::desc(amplitude_mean)) %>%
                       dplyr::mutate(amplitude_mean_rank = row_number()) %>%
                       as.data.frame() %>%
                       magrittr::set_rownames(.$gene)
  if (verbose) {
    cli::cli_alert_success("Posterior summarization complete.")
    cli::cli_alert_info("Mean amplitude of control data: ", round(mean(amplitude_summary$amplitude_mean), 4))
  }
  return(amplitude_summary)
}
