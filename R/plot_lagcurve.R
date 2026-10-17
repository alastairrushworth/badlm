#' Plot the estimated lag curve
#'
#' Draws the posterior median of the lag curve together with its pointwise 95\%
#' credible interval.
#'
#' @param dlm_object A fitted model returned by [badlm()].
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' truth  <- dnorm(0:20, mean = 6, sd = 3) / 10
#' expose <- as.numeric(arima.sim(list(ar = 0.5), n = 300))
#' resp   <- lag_matrix(expose, p = 20) %*% truth + rnorm(280, sd = 0.01)
#' fit    <- badlm(expose, resp, k = 10, nlag = 20, samples = 1000, verbose = FALSE)
#' plot_lagcurve(fit)
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab ggtitle
#' @export
plot_lagcurve <- function(dlm_object) {
  cf <- dlm_object$coefs
  lagdf <- tibble::tibble(
    lag      = rep(cf$lag, 3),
    quantile = rep(c("2.5%", "median", "97.5%"), each = nrow(cf)),
    coef     = c(cf$lower, cf$beta, cf$upper)
  )
  ggplot(lagdf, aes(x = lag, y = coef, group = quantile, colour = quantile)) +
    geom_line() +
    xlab("Lag number") +
    ylab("Lag coefficient") +
    ggtitle("Posterior median and credible interval for lag function")
}
