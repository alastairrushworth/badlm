#' badlm: Bayesian adaptive distributed lag models
#'
#' Fits the adaptive distributed lag model of Rushworth (2018). The lag curve
#' is represented by a B-spline basis whose coefficients receive a
#' first-difference penalty with lag-varying precision; the log precisions are
#' themselves smoothed by a random walk prior, so that the amount of smoothing
#' adapts to the shape of the lag curve. See [badlm()] to fit a model and
#' [plot_lagcurve()] to display the estimated lag curve.
#'
#' @references Rushworth, A. (2018). Bayesian adaptive distributed lag models.
#'   arXiv:1801.06670. \url{https://arxiv.org/abs/1801.06670}
#' @keywords internal
"_PACKAGE"

# Column names used inside ggplot2::aes() in plot_lagcurve()
utils::globalVariables(c("lag", "coef", "quantile"))
