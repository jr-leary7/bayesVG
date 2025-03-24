#' A \code{ggplot2} theme for \code{bayesVG}.
#'
#' @name theme_bayesVG
#' @author Jack R. Leary
#' @description A publication-ready theme for creating gene expression scatterplots, embedding plots, etc.
#' @param base.size The base font size. Defaults to 12.
#' @param base.lwd The base linewidth. Defaults to 0.75.
#' @param base.family The font family to be used throughout. Defaults to "sans".
#' @param umap (Optional) If set to TRUE, removes axis text and ticks for a cleaner look. Defaults to FALSE.
#' @param spatial (Optional) If set to true, provides a box around the spatial representation of the data and cleans up the axes. Defaults to FALSE.
#' @importFrom cli cli_abort
#' @importFrom ggplot2 theme_classic theme element_rect element_line element_blank
#' @return A \code{ggplot2} theme.
#' @export

theme_bayesVG <- function(base.size = 12,
                          base.lwd = 0.75,
                          base.family = "sans",
                          umap = FALSE,
                          spatial = FALSE) {
  # check inputs
  if (umap && spatial) { cli::cli_abort("Only one of {.field umap} and {.field spatial} can be specified at once.") }
  # generate theme
  theme_bayesVG <- ggplot2::theme_classic(base_size = base.size,
                                          base_family = base.family,
                                          base_line_size = base.lwd,
                                          base_rect_size = base.lwd) +
                   ggplot2::theme(strip.clip = "off",
                                  strip.background = ggplot2::element_rect(linewidth = base.lwd),
                                  axis.line = ggplot2::element_line(lineend = "square"))
  if (umap) {
    theme_bayesVG <- theme_bayesVG +
                     ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                    axis.text = ggplot2::element_blank())
  } else if (spatial) {
    theme_bayesVG <- theme_bayesVG +
                     ggplot2::theme(axis.line = ggplot2::element_blank(),
                                    axis.ticks = ggplot2::element_blank(),
                                    axis.text = ggplot2::element_blank(),
                                    panel.border = ggplot2::element_rect(colour = "black",
                                                                         fill = NA,
                                                                         linewidth = 2 * base.lwd))
  }
  return(theme_bayesVG)
}
