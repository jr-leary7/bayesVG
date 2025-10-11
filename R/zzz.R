.onLoad <- function(...) {
  cmdstan_version <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
  if (is.null(cmdstan_version)) {
    cli::cli_abort("No CmdStan installation found. Run cmdstanr::install_cmdstan() to install.")
  }
}

utils::globalVariables(c(
  "gene", 
  "amplitude_mean_rank",
  "amplitude_mean", 
  "i",
  "r", 
  "p_value",
  ".", 
  "spot", 
  "gene_expression", 
  "x", 
  "y", 
  "g", 
  "amplitude_q5", 
  "amplitude_q95", 
  "variable",
  "gene_id", 
  "amplitude_sd", 
  "subject", 
  "quintile", 
  "cell", 
  "b_Intercept", 
  "intercept", 
  "gene_re", 
  "subject_re", 
  "mu", 
  "b_shape_Intercept",
  "theta_re", 
  "theta", 
  "sigma2", 
  "dispersion", 
  "dispersion_mean", 
  "hvg_status", 
  "mu_mean", 
  "hvg_label", 
  "module_score", 
  "dim1", 
  "dim2",
  "amplitude_mean_rank", 
  "svg_status", 
  "mu_naive", 
  "svg_label", 
  "meta_vec", 
  "gene_expr", 
  "b", 
  "hex"
))
