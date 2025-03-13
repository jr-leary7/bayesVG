.onLoad <- function(...) {
  cmdstan_version <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
  if (is.null(cmdstan_version)) {
    cli::cli_abort("No CmdStan installation found. Run cmdstanr::install_cmdstan() to install.")
  }
}
