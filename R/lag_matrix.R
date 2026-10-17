#' Matrix of lagged values of a time series
#'
#' Builds the design matrix of a distributed lag model. Row \eqn{t} contains
#' \eqn{x_t, x_{t-1}, \ldots, x_{t-p}} for \eqn{t = p + 1, \ldots, n}; the
#' first \eqn{p} time points, for which the full set of lags is unavailable,
#' are dropped.
#'
#' @param x Numeric vector, the time series to lag.
#' @param p Integer, the maximum lag.
#' @param start.at.zero Logical. If `TRUE` (the default) the first column is
#'   lag 0, that is \eqn{x_t} itself, and the matrix has \eqn{p + 1} columns.
#'   Otherwise the columns start at lag 1 and the matrix has \eqn{p} columns.
#' @return A numeric matrix with `length(x) - p` rows.
#' @examples
#' lag_matrix(1:6, p = 2)
#' lag_matrix(1:6, p = 2, start.at.zero = FALSE)
#' @export
lag_matrix <- function(x, p, start.at.zero = TRUE) {
  x <- as.numeric(x)
  p <- as.integer(p)
  if (length(p) != 1L || is.na(p) || p < 1L) {
    stop("`p` must be a single positive integer")
  }
  if (p >= length(x)) {
    stop("`p` must be smaller than `length(x)`")
  }
  m <- stats::embed(x, p + 1L)
  if (!start.at.zero) m <- m[, -1L, drop = FALSE]
  m
}
