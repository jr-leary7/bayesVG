#' Identify highly variable genes in a Bayesian manner.
#'
#' @name findVariableFeaturesBayes
#' @author Jack R. Leary
#' @description This function implements HVG estimation using Bayesian variational inference to approximate the posterior distribution of the mean and overdispersion of each gene.
#' @param sc.obj An object of class \code{Seurat} or \code{SingleCellExperiment}. Defaults to NULL.
#' @param subject.id A string specifying the metadata column in \code{sc.obj} that contains subject IDs. Defaults to NULL.
#' @param n.cells.subsample An integer specifying the number of cells per-gene (and per-subject if \code{subject.id} is not NULL) to subsample to when performing estimation. Defaults to 500.
#' @param iter.per.chain An integer specifying the number of iterations per chain. Defaults to 2000.
#' @param warmup.per.chain An integer specifying the number of warmup (burn-in) iterations per chain. Defaults to 250.
#' @param n.chains An integer specifying the number of chains used when fitting the model via sampling instead of VI or the Laplace approximation. Defaults to 1.
#' @param n.cores.chain An integer specifying the number of cores to be used when fitting the Bayesian hierarchical model. Defaults to 1.
#' @param n.cores.per.chain An integer specifying the number of cores to be used within each chain when fitting the Bayesian hierarchical model. Defaults to 4.
#' @param model.priors A vector containing priors to be used in model fitting. If left NULL, intelligent priors will be set internally. See \code{\link[brms]{set_prior}} for details. Defaults to NULL.
#' @param algorithm A string specifying the variational inference or sampling algorithm to be used. Must be one of "meanfield", "fullrank", "laplace", "pathfinder" (all approximate methods) or "sampling" for MCMC via NUTS. Note that MCMC will be slower than using a VI algorithm or the Laplace approximation. Defaults to "meanfield".
#' @param opencl.params A two-element double vector specifying the platform and device IDs of the OpenCL GPU device. Most users should specify \code{c(0, 0)}. See \code{\link[brms]{opencl}} for more details. Defaults to NULL.
#' @param random.seed A double specifying the random seed to be used when fitting and sampling from the model. Defaults to 312.
#' @param verbose A Boolean specifying whether or not verbose model output should be printed to the console. Defaults to TRUE.
#' @param save.model A Boolean specifying whether or not the fitted model generated by \code{\link[brms]{brm}} should be saved to the unstructured metadata of \code{sc.obj}. Defaults to FALSE.
#' @details
#' \itemize{
#' \item Throughout the package, we make an important distinction between overdispersion (the parameter \eqn{\theta} of the Negative-binomial distribution) and dispersion, which is estimated as \eqn{d = \frac{\sigma^2}{\mu}}.
#' \item Our method makes use of \code{cmdstanr} to fit the model rather than \code{rstan}, as the former is generally much faster. For details, see \code{\link[cmdstanr]{cmdstan_model}}. This of course necessitates first running \code{\link[cmdstanr]{install_cmdstan}} if you haven't already. If errors occur, make sure to check that your toolchain is set up correctly by running \code{\link[cmdstanr]{check_cmdstan_toolchain}}.
#' \item When using sampling instead of an approximation method, increasing \code{n.chains} will increase the model's performance at the cost of extra computational resource usage. If possible, set \code{n.cores} equal to \code{n.chains} for optimal processing speed.
#' \item While we have implemented GPU acceleration via OpenCL through the argument \code{opencl.params}, OpenCL acceleration is not supported on every machine. For example, Apple M-series chips do not support double-precision floating-points, which are necessary for Stan to compile. For more information, see \href{https://discourse.mc-stan.org/t/gpus-on-mac-osx-apple-m1/23375/5}{this Stan forums thread}. In order to correctly specify the OpenCL platform and device IDs, use the \code{clinfo} command line utility.
#' \item The user can specify which variational inference (VI) or sampling algorithm to use to fit the model via the argument \code{algorithm}. For further details, see \href{https://www.jmlr.org/papers/volume18/16-107/16-107.pdf}{this paper} comparing the meanfield and fullrank algorithms, and \href{https://doi.org/10.48550/arXiv.2108.03782}{this preprint} that introduced the Pathfinder algorithm. For a primer on automatic differentiation variational inference (ADVI), see \href{https://doi.org/10.48550/arXiv.1506.03431}{this preprint}. Lastly, \href{https://discourse.mc-stan.org/t/issues-with-differences-between-mcmc-and-pathfinder-results-how-to-make-pathfinder-or-something-else-more-accurate/35992}{this Stan forums thread} lays out some pratical differences between the algorithms.
#' \item If \code{save.model} is set to TRUE, the final model fit will be saved to the appropriate unstructured metadata slot of \code{sc.obj}. This allows the user to inspect the final fit and perform posterior predictive checks, but the model object takes up a lot of space. As such, it is recommended to remove it from \code{sc.obj} by setting the appropriate slot to NULL before saving it to disk.
#' }
#' @import cmdstanr
#' @import magrittr
#' @importFrom cli cli_abort cli_alert_warning cli_alert_success
#' @importFrom parallelly availableCores
#' @importFrom SingleCellExperiment colData logcounts
#' @importFrom SummarizedExperiment rowData
#' @importFrom BiocGenerics counts
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom dplyr mutate select with_groups ntile slice_sample n filter summarise arrange desc left_join inner_join pull row_number
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches all_of everything
#' @importFrom stats quantile var
#' @importFrom withr with_output_sink
#' @importFrom brms set_prior brm bf negbinomial
#' @importFrom posterior as_draws_df
#' @importFrom S4Vectors DataFrame
#' @importFrom methods slot
#' @importFrom purrr map
#' @return Depending on the input, either an object of class \code{Seurat} or \code{SingleCellExperiment} with gene-level statistics added to the appropriate metadata slot.
#' @seealso \code{\link[Seurat]{FindVariableFeatures}}
#' @seealso \code{\link[scran]{modelGeneVar}}
#' @seealso \code{\link[brms]{brm}}
#' @export
#' @examples
#' data(seu_pbmc)
#' seu_pbmc <- findVariableFeaturesBayes(seu_pbmc,
#'                                       n.cells.subsample = 500L,
#'                                       algorithm = "meanfield",
#'                                       n.cores.per.chain = 1L,
#'                                       save.model = TRUE)

findVariableFeaturesBayes <- function(sc.obj = NULL,
                                      subject.id = NULL,
                                      n.cells.subsample = 500L,
                                      iter.per.chain = 2000L,
                                      warmup.per.chain = 250L,
                                      n.chains = 1L,
                                      n.cores.chain = 1L,
                                      n.cores.per.chain = 4L,
                                      model.priors = NULL,
                                      algorithm = "meanfield",
                                      opencl.params = NULL,
                                      random.seed = 312,
                                      verbose = TRUE,
                                      save.model = FALSE) {
  # check & parse inputs
  if (is.null(sc.obj)) { cli::cli_abort("Please provide all inputs to findVariableFeaturesBayes().") }
  n_cores_total <- n.cores.chain * n.cores.per.chain
  if (n_cores_total > unname(parallelly::availableCores())) { cli::cli_abort("The total number of requested cores is greater than the number of available cores on your machine.") }
  if (n.cores.chain != n.chains) { cli::cli_alert_warning("In general, the number of cores should equal the number of chains for optimal performance.") }
  algorithm <- tolower(algorithm)
  if (!algorithm %in% c("meanfield", "fullrank", "pathfinder", "laplace", "sampling")) { cli::cli_abort("Please provide a valid sampling or approximation algorithm.") }
  if (algorithm == "sampling" && n.chains == 1L) { cli::cli_alert_warning("It is recommended to use multiple chains when utilizing MCMC sampling.") }
  if (!is.null(opencl.params) && (!is.double(opencl.params) || !length(opencl.params) == 2)) { cli::cli_abort("Argument {.field opencl.params} must be a double vector of length 2 if non-NULL.") }
  if (is.null(opencl.params)) {
    opencl_IDs <- NULL
  } else {
    opencl_IDs <- opencl.params
  }
  if (is.null(model.priors)) {
    model.priors <- c(brms::set_prior("normal(0, 2)", class = "Intercept", dpar = "mu"),
                      brms::set_prior("student_t(3, 0, 2)", class = "sd", dpar = "mu"),
                      brms::set_prior("normal(0, 2)", class = "Intercept", dpar = "shape"),
                      brms::set_prior("student_t(3, 0, 2)", class = "sd", dpar = "shape"))
  }
  # start time tracking
  time_start <- Sys.time()
  # extract (sparse) counts matrix
  if (inherits(sc.obj, "SingleCellExperiment")) {
    if (!is.null(subject.id)) {
      subject_vec <- SingleCellExperiment::colData(sc.obj)[[subject.id]]
    }
    expr_mat <- BiocGenerics::counts(sc.obj)
  } else if (inherits(sc.obj, "Seurat")) {
    if (!is.null(subject.id)) {
      subject_vec <- sc.obj@meta.data[[subject.id]]
    }
    expr_mat <- Seurat::GetAssayData(sc.obj,
                                     layer = "counts",
                                     assay = Seurat::DefaultAssay(sc.obj))
  }
  # generate offset from counts matrix 
  seq_depth <- log(Matrix::colSums(expr_mat))
  # convert counts matrix to long data.frame for modeling
  expr_df <- as.data.frame(expr_mat) %>%
             dplyr::mutate(gene = rownames(.), .before = 1)
  sampled_cells_per_quintile <- round(n.cells.subsample / 5)
  if (!is.null(subject.id)) {
    expr_df <- tidyr::pivot_longer(expr_df,
                                   cols = !gene,
                                   names_to = "cell",
                                   values_to = "gene_expression") %>%
               dplyr::mutate(subject = rep(subject_vec, times = nrow(.) / length(subject_vec)), .before = 1) %>%
               dplyr::with_groups(c(gene, subject),
                                  dplyr::mutate,
                                  quintile = dplyr::ntile(gene_expression, 5)) %>%
               dplyr::with_groups(c(gene, subject, quintile),
                                  dplyr::slice_sample,
                                  n = sampled_cells_per_quintile) %>%
               dplyr::mutate(gene = factor(gene, levels = unique(gene)),
                             subject = factor(subject, levels = unique(subject)))
  } else {
    expr_df <- tidyr::pivot_longer(expr_df,
                                   cols = !gene,
                                   names_to = "cell",
                                   values_to = "gene_expression") %>%
               dplyr::with_groups(gene,
                                  dplyr::mutate,
                                  quintile = dplyr::ntile(gene_expression, 5)) %>%
               dplyr::with_groups(c(gene, quintile),
                                  dplyr::slice_sample,
                                  n = sampled_cells_per_quintile) %>%
               dplyr::mutate(gene = factor(gene, levels = unique(gene)))
  }
  # convert from tibble to data.frame & convert gene expression to integer to save space
  expr_df <- as.data.frame(expr_df) %>%
             dplyr::mutate(gene_expression = as.integer(gene_expression)) %>%
             dplyr::select(-c(cell, quintile))
  # create model formula
  if (!is.null(subject.id)) {
    model_formula <- brms::bf(gene_expression ~ 1 + (1 | gene) + (1 | subject),
                              shape ~ 1 + (1 | gene))
  } else {
    model_formula <- brms::bf(gene_expression ~ 1 + (1 | gene),
                              shape ~ 1 + (1 | gene))
  }
  # fit negative-binomial hierarchical bayesian model
  if (verbose) {
    brms_fit <- brms::brm(model_formula,
                          prior = model.priors,
                          data = expr_df,
                          family = brms::negbinomial(link = "log", link_shape = "log"),
                          chains = n.chains,
                          iter = iter.per.chain,
                          warmup = warmup.per.chain,
                          cores = n.cores.chain,
                          threads = n.cores.per.chain,
                          opencl = opencl_IDs,
                          normalize = FALSE,
                          silent = 2,
                          backend = "cmdstanr",
                          algorithm = algorithm,
                          stan_model_args = list(stanc_options = list("O1")),
                          seed = random.seed)
  } else {
    withr::with_output_sink(tempfile(), {
      brms_fit <- brms::brm(model_formula,
                            prior = model.priors,
                            data = expr_df,
                            family = brms::negbinomial(link = "log", link_shape = "log"),
                            chains = n.chains,
                            iter = iter.per.chain,
                            warmup = warmup.per.chain,
                            cores = n.cores.chain,
                            threads = n.cores.per.chain,
                            opencl = opencl_IDs,
                            normalize = FALSE,
                            silent = 2,
                            backend = "cmdstanr",
                            algorithm = algorithm,
                            stan_model_args = list(stanc_options = list("O1")),
                            seed = random.seed)
    })
  }
  # draw samples from approximate posterior
  posterior_samples <- as.data.frame(posterior::as_draws_df(brms_fit))
  # estimate posterior gene means
  mu_intercept <- dplyr::pull(posterior_samples, b_Intercept)
  gene_random_effects <- dplyr::select(posterior_samples, tidyselect::matches("r_gene\\[.*Intercept")) %>%
                         dplyr::mutate(intercept = mu_intercept,
                                       sample = dplyr::row_number(),
                                       .before = 1)
  if (!is.null(subject.id)) {
    subject_random_effects <- dplyr::select(posterior_samples, tidyselect::matches("r_subject\\[.*Intercept")) %>%
                              dplyr::mutate(sample = dplyr::row_number(), .before = 1)
    subject_samples_long <- tidyr::pivot_longer(subject_random_effects,
                                                cols = !sample,
                                                names_to = "subject",
                                                values_to = "subject_re") %>%
                            as.data.frame() %>%
                            dplyr::mutate(subject = gsub(",Intercept\\]", "", gsub("r_subject\\[", "", subject)))
    gene_samples_long <- tidyr::pivot_longer(gene_random_effects,
                                             cols = !c(intercept, sample),
                                             names_to = "gene",
                                             values_to = "gene_re") %>%
                         as.data.frame() %>%
                         dplyr::mutate(gene = gsub(",Intercept\\]", "", gsub("r_gene\\[", "", gene)))
    mu_samples_long <- dplyr::inner_join(gene_samples_long,
                                         subject_samples_long,
                                         by = "sample",
                                         relationship = "many-to-many") %>%
                       dplyr::mutate(mu = exp(intercept + gene_re + subject_re)) %>%
                       dplyr::select(-intercept)
    mu_summary <- dplyr::with_groups(mu_samples_long,
                                     c(gene, subject),
                                     dplyr::summarise,
                                     mu_mean = mean(mu),
                                     mu_var = stats::var(mu),
                                     mu_ci_ll = stats::quantile(mu, 0.025),
                                     mu_ci_ul = stats::quantile(mu, 0.975))
  } else {
    mu_samples_long <- tidyr::pivot_longer(gene_random_effects,
                                           cols = !c(intercept, sample),
                                           names_to = "gene",
                                           values_to = "gene_re") %>%
                       as.data.frame() %>%
                       dplyr::mutate(gene = gsub(",Intercept\\]", "", gsub("r_gene\\[", "", gene)),
                                     mu = exp(intercept + gene_re)) %>%
                       dplyr::select(-intercept)
    mu_summary <- dplyr::with_groups(mu_samples_long,
                                     gene,
                                     dplyr::summarise,
                                     mu_mean = mean(mu),
                                     mu_var = stats::var(mu),
                                     mu_ci_ll = stats::quantile(mu, 0.025),
                                     mu_ci_ul = stats::quantile(mu, 0.975))
  }
  # estimate posterior gene overdispersions
  theta_intercept <- dplyr::pull(posterior_samples, b_shape_Intercept)
  theta_random_effects <- dplyr::select(posterior_samples, tidyselect::matches("r_gene__shape\\[.*Intercept")) %>%
                          dplyr::mutate(intercept = theta_intercept,
                                        sample = dplyr::row_number(),
                                        .before = 1)
  theta_samples_long <- tidyr::pivot_longer(theta_random_effects,
                                            cols = !c(intercept, sample),
                                            names_to = "gene",
                                            values_to = "theta_re") %>%
                        as.data.frame() %>%
                        dplyr::mutate(gene = gsub(",Intercept\\]", "", gsub("r_gene__shape\\[", "", gene)),
                                      theta = exp(intercept + theta_re)) %>%
                        dplyr::select(-intercept)
  theta_summary <- dplyr::with_groups(theta_samples_long,
                                      gene,
                                      dplyr::summarise,
                                      theta_mean = mean(theta),
                                      theta_var = stats::var(theta),
                                      theta_ci_ll = stats::quantile(theta, 0.025),
                                      theta_ci_ul = stats::quantile(theta, 0.975))
  # estimate posterior gene variances & dispersions
  if (!is.null(subject.id)) {
    sigma2_samples_long <- dplyr::inner_join(mu_samples_long,
                                             theta_samples_long,
                                             by = c("gene", "sample"),
                                             relationship = "many-to-many") %>%
                           dplyr::mutate(sigma2 = mu * (1 + mu / theta),
                                         dispersion = sigma2 / mu)
    sigma2_summary <- dplyr::with_groups(sigma2_samples_long,
                                         c(gene, subject),
                                         dplyr::summarise,
                                         sigma2_mean = mean(sigma2),
                                         sigma2_var = stats::var(sigma2),
                                         sigma2_ci_ll = stats::quantile(sigma2, 0.025),
                                         sigma2_ci_ul = stats::quantile(sigma2, 0.975),
                                         dispersion_mean = mean(dispersion),
                                         dispersion_var = var(dispersion),
                                         dispersion_ci_ll = stats::quantile(dispersion, 0.025),
                                         dispersion_ci_ul = stats::quantile(dispersion, 0.975))
    gene_summary <- dplyr::inner_join(mu_summary,
                                      theta_summary,
                                      by = "gene") %>%
                    dplyr::inner_join(sigma2_summary, by = c("gene", "subject"))
  } else {
    sigma2_samples_long <- dplyr::inner_join(mu_samples_long,
                                             theta_samples_long,
                                             by = c("gene", "sample")) %>%
                           dplyr::mutate(sigma2 = mu * (1 + mu / theta),
                                         dispersion = sigma2 / mu)
    sigma2_summary <- dplyr::with_groups(sigma2_samples_long,
                                         gene,
                                         dplyr::summarise,
                                         sigma2_mean = mean(sigma2),
                                         sigma2_var = stats::var(sigma2),
                                         sigma2_ci_ll = stats::quantile(sigma2, 0.025),
                                         sigma2_ci_ul = stats::quantile(sigma2, 0.975),
                                         dispersion_mean = mean(dispersion),
                                         dispersion_var = var(dispersion),
                                         dispersion_ci_ll = stats::quantile(dispersion, 0.025),
                                         dispersion_ci_ul = stats::quantile(dispersion, 0.975))
    gene_summary <- dplyr::inner_join(mu_summary,
                                      theta_summary,
                                      by = "gene") %>%
                    dplyr::inner_join(sigma2_summary, by = "gene") %>%
                    magrittr::set_rownames(.$gene)
  }
  if (verbose) {
    cli::cli_alert_success("Posterior summarization complete.")
  }
  # add gene-level estimates to object metadata
  if (is.null(subject.id)) {
    if (inherits(sc.obj, "SingleCellExperiment")) {
      gene_summary_s4 <- SummarizedExperiment::rowData(sc.obj) %>%
                         as.data.frame() %>%
                         dplyr::mutate(gene = rownames(.), .before = 1) %>%
                         dplyr::left_join(gene_summary, by = "gene") %>%
                         S4Vectors::DataFrame()
      SummarizedExperiment::rowData(sc.obj) <- gene_summary_s4
    } else if (inherits(sc.obj, "Seurat")) {
      version_check <- try({
        methods::slot(sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]], name = "meta.data")
      }, silent = TRUE)
      if (inherits(version_check, "try-error")) {
        orig_metadata <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features
      } else {
        orig_metadata <- sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data
      }
      if (ncol(orig_metadata) > 0) {
        new_metadata <- dplyr::mutate(orig_metadata,
                                      gene = rownames(sc.obj),
                                      .before = 1) %>%
                        dplyr::left_join(gene_summary, by = "gene")
        if (inherits(version_check, "try-error")) {
          sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features <- new_metadata
        } else {
          sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- new_metadata
        }
      } else {
        gene_summary <- gene_summary[rownames(sc.obj), ]
        if (inherits(version_check, "try-error")) {
          sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.features <- gene_summary
        } else {
          sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@meta.data <- gene_summary
        }
      }
    }
  } else {
    gene_summary_list <- split(gene_summary, gene_summary$subject)
    gene_summary_list <- purrr::map(gene_summary_list, \(x) {
      rownames(x) <- x$gene
      x <- x[rownames(sc.obj), ]
      return(x)
    })
    if (inherits(sc.obj, "SingleCellExperiment")) {
      sc.obj@metadata$gene_stats_bayes <- gene_summary_list
    } else if (inherits(sc.obj, "Seurat")) {
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@misc$gene_stats_bayes <- gene_summary_list
    }
  }
  # optionally save model fit to object's unstructured metadata
  if (save.model) {
    if (inherits(sc.obj, "SingleCellExperiment")) {
      sc.obj@metadata$model_fit <- brms_fit
    } else if (inherits(sc.obj, "Seurat")) {
      sc.obj@assays[[Seurat::DefaultAssay(sc.obj)]]@misc$model_fit <- brms_fit
    }
  }
  # finish time tracking
  time_diff <- Sys.time() - time_start
  time_units <- ifelse(attributes(time_diff)$units == "secs",
                       "seconds",
                       ifelse(attributes(time_diff)$units == "mins",
                              "minutes",
                              "hours"))
  if (verbose) {
    time_message <- paste0("{.pkg bayesVG} modeling of ",
                           nrow(sc.obj),
                           " genes completed in ",
                           as.numeric(round(time_diff, 3)),
                           " ",
                           time_units,
                           ".")
    cli::cli_alert_success(time_message)
  }
  return(sc.obj)
}
