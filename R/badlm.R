#' Fit a Bayesian adaptive distributed lag model
#'
#' Estimates the distributed lag relationship between an exposure series `x`
#' and a response series `y`, following Rushworth (2018). The lag curve over
#' lags \eqn{0, \ldots, p} is represented by `k` cubic B-spline basis
#' functions. The basis coefficients receive a first-difference penalty whose
#' precision varies along the lag axis, and the log precisions are smoothed by
#' a random walk prior, so that the amount of smoothing adapts to the shape of
#' the lag curve. The curve is softly constrained to be zero at the maximum lag
#' through a single pseudo-observation. Inference is by Markov chain Monte
#' Carlo.
#'
#' @details The sampler alternates Gibbs updates of the basis coefficients, the
#'   error variance and the random walk variance with one-at-a-time random-walk
#'   Metropolis updates of the log precisions. Metropolis proposal scales are
#'   tuned during the burn-in period only, so the retained draws form a valid
#'   Markov chain. Credible intervals are computed pointwise from posterior
#'   draws of the lag curve itself.
#'
#' @param x Numeric vector, the exposure or predictor series.
#' @param y Numeric vector, the response series. Either the same length as
#'   `x`, in which case the first `nlag` values are dropped, or of length
#'   `length(x) - nlag`.
#' @param k Integer, the number of B-spline basis functions (at least 4).
#' @param nlag Integer, the maximum lag \eqn{p}.
#' @param samples Total number of MCMC iterations.
#' @param burnin Number of initial iterations to discard. Defaults to half of
#'   `samples`.
#' @param thin Thinning interval for the retained draws.
#' @param verbose Logical; print a progress bar?
#'
#' @return A list with components
#' \describe{
#'   \item{`coefs`}{A tibble with one row per lag holding the posterior median
#'     (`beta`) and the 2.5\% and 97.5\% posterior quantiles (`lower`, `upper`)
#'     of the lag curve.}
#'   \item{`lag_samples`}{Matrix of posterior draws of the lag curve, one row
#'     per retained iteration and one column per lag.}
#'   \item{`ed`}{Posterior draws of the effective degrees of freedom of the
#'     fitted lag curve.}
#'   \item{`model`}{The raw sampler output: draws of the basis coefficients
#'     (`bet_store`), the lag-varying precisions (`lam_store`), the error
#'     variance (`sig_store`) and the random walk variance (`zeta_store`),
#'     together with the Metropolis acceptance rates.}
#'   \item{`basis`}{The B-spline basis matrix, with one row per lag.}
#' }
#' @references Rushworth, A. (2018). Bayesian adaptive distributed lag models.
#'   arXiv:1801.06670. \url{https://arxiv.org/abs/1801.06670}
#' @seealso [plot_lagcurve()], [lag_matrix()]
#' @examples
#' set.seed(1)
#' truth  <- dnorm(0:20, mean = 6, sd = 3) / 10
#' expose <- as.numeric(arima.sim(list(ar = 0.5), n = 300))
#' resp   <- lag_matrix(expose, p = 20) %*% truth + rnorm(280, sd = 0.01)
#' fit    <- badlm(expose, resp, k = 10, nlag = 20, samples = 1000, verbose = FALSE)
#' fit$coefs
#' @export
badlm <- function(x, y, k, nlag, samples = 10000, burnin = floor(samples / 2),
                  thin = 1, verbose = TRUE) {
  x    <- as.numeric(x)
  y    <- as.numeric(y)
  nlag <- as.integer(nlag)
  k    <- as.integer(k)
  if (anyNA(x) || anyNA(y)) stop("missing values in `x` or `y` are not supported")
  if (nlag < 1L || nlag >= length(x)) stop("`nlag` must be between 1 and `length(x)` - 1")
  if (k < 4L) stop("`k` must be at least 4")
  if (length(y) == length(x)) y <- y[-seq_len(nlag)]
  if (length(y) != length(x) - nlag) {
    stop("`y` must have length `length(x)` or `length(x) - nlag`")
  }

  lag_mat <- lag_matrix(x, p = nlag)
  bbase   <- b_spline_basis(0:nlag, nseg = k - 3L, range.variables = c(0, nlag))
  # One pseudo-observation with zero response pulls the curve towards zero at lag `nlag`.
  X_mat   <- rbind(lag_mat %*% bbase, bbase[nlag + 1L, ])
  y_aug   <- c(y, 0)

  model_out <- gaussian_badlm(X = X_mat, y = y_aug, samples = samples,
                              burnin = burnin, thin = thin, verbose = verbose)

  # posterior draws of the lag curve itself, then pointwise quantiles
  lag_samples <- model_out$bet_store %*% t(bbase)
  qs <- apply(lag_samples, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  coefs <- tibble::tibble(lag   = 0:nlag,
                          beta  = unname(qs[2, ]),
                          upper = unname(qs[3, ]),
                          lower = unname(qs[1, ]))

  list(coefs       = coefs,
       lag_samples = lag_samples,
       ed          = model_out$ed_store,
       model       = model_out,
       basis       = bbase)
}
