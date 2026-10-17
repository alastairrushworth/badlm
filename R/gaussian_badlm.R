# MCMC sampler for the Gaussian adaptive distributed lag model.
#
# In the notation of Rushworth (2018, Section 2) the model is
#
#   y | b, sigma2              ~ N(X b, sigma2 I)
#   b | lambda                 ~ N(0, P^-),      P = D' diag(lambda) D   (intrinsic, rank K - 1)
#   tau = log(lambda) | zeta2  ~ N(0, zeta2 K^-)                         (intrinsic random walk, rank K - 2)
#   sigma2 ~ IG(a_sigma, b_sigma),  zeta2 ~ IG(a_zeta, b_zeta)
#
# where D is the (K - 1) x K first-difference matrix. b, sigma2 and zeta2 have
# closed-form full conditionals and are updated by Gibbs steps; each tau_j is
# updated by a random-walk Metropolis step whose scale is tuned during burn-in.
#
# P is singular, so the normalising constant of pi(b | lambda) involves the
# generalised determinant of P, the product of its non-zero eigenvalues. For the
# chain graph the matrix-tree theorem gives this in closed form,
#
#   pdet(P) = K * prod(lambda_k),
#
# so that  log pi(b | lambda) = 0.5 * sum(tau) - 0.5 * sum(lambda * diff(b)^2) + const.
# Earlier versions of this function evaluated determinant(P) numerically, which
# for a singular matrix returns the logarithm of rounding error.

gaussian_badlm <- function(X, y, samples, burnin = floor(samples / 2), thin = 1,
                           prior_sigma2 = c(0.001, 0.001), prior_zeta2 = c(1, 0.5),
                           tau_max = 40, tune_every = 200, verbose = TRUE) {
  samples <- as.integer(samples)
  burnin  <- as.integer(burnin)
  thin    <- as.integer(thin)
  if (samples < 1L) stop("`samples` must be positive")
  if (burnin < 0L || burnin >= samples) stop("`burnin` must lie between 0 and `samples` - 1")
  if (thin < 1L) stop("`thin` must be a positive integer")
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (nrow(X) != length(y)) stop("`X` and `y` have incompatible dimensions")
  K <- ncol(X)
  n <- nrow(X)
  if (K < 3L) stop("at least three basis functions are required")

  XTX <- crossprod(X)
  Xty <- as.numeric(crossprod(X, y))
  D   <- diff(diag(K))                      # (K - 1) x K first differences

  n_save     <- (samples - burnin) %/% thin
  bet_store  <- matrix(0, n_save, K)
  lam_store  <- matrix(0, n_save, K - 1L)
  sig_store  <- numeric(n_save)
  zeta_store <- numeric(n_save)
  ed_store   <- numeric(n_save)

  # initial values: start from a heavily smoothed curve
  beta     <- numeric(K)
  sigma2   <- stats::var(y)
  zeta2    <- 0.01
  tau      <- stats::rnorm(K - 1L, mean = 8, sd = 1)
  sd_prop  <- rep(1, K - 1L)
  acc_tune <- numeric(K - 1L)     # acceptances since the last tuning step
  acc_post <- numeric(K - 1L)     # acceptances after burn-in

  if (verbose) {
    cat("Collecting", samples, "samples\n")
    pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  }
  saved <- 0L
  for (i in seq_len(samples)) {

    ## tau_j | b, zeta2: one-at-a-time random-walk Metropolis
    d2 <- diff(beta)^2
    for (j in seq_len(K - 1L)) {
      tau_new <- stats::rnorm(1, tau[j], sd_prop[j])
      if (tau_new > tau_max) next
      tau_prop    <- tau
      tau_prop[j] <- tau_new
      log_ratio <- 0.5 * (tau_new - tau[j]) -
        0.5 * (exp(tau_new) - exp(tau[j])) * d2[j] -
        (sum(diff(tau_prop)^2) - sum(diff(tau)^2)) / (2 * zeta2)
      if (log(stats::runif(1)) < log_ratio) {
        tau[j]      <- tau_new
        acc_tune[j] <- acc_tune[j] + 1
        if (i > burnin) acc_post[j] <- acc_post[j] + 1
      }
    }
    if (i <= burnin && i %% tune_every == 0L) {
      rate       <- acc_tune / tune_every
      sd_prop    <- sd_prop * ifelse(rate > 0.4, 2, ifelse(rate < 0.3, 0.5, 1))
      acc_tune[] <- 0
    }

    ## b | tau, sigma2 ~ N(prec^{-1} m, prec^{-1}), prec = P + X'X / sigma2, m = X'y / sigma2
    lambda <- exp(tau)
    P      <- crossprod(sqrt(lambda) * D)
    R      <- chol(P + XTX / sigma2)
    mu     <- backsolve(R, backsolve(R, Xty / sigma2, transpose = TRUE))
    beta   <- as.numeric(mu + backsolve(R, stats::rnorm(K)))

    ## sigma2 | b
    rss    <- sum((y - X %*% beta)^2)
    sigma2 <- 1 / stats::rgamma(1, shape = prior_sigma2[1] + n / 2,
                                rate = prior_sigma2[2] + rss / 2)

    ## zeta2 | tau   (the random walk prior on tau has rank K - 2)
    zeta2  <- 1 / stats::rgamma(1, shape = prior_zeta2[1] + (K - 2) / 2,
                                rate = prior_zeta2[2] + sum(diff(tau)^2) / 2)

    if (i > burnin && (i - burnin) %% thin == 0L) {
      saved <- saved + 1L
      bet_store[saved, ] <- beta
      lam_store[saved, ] <- lambda
      sig_store[saved]   <- sigma2
      zeta_store[saved]  <- zeta2
      # effective degrees of freedom, tr[X (X'X + sigma2 P)^{-1} X']
      ed_store[saved]    <- sum(diag(solve(XTX + sigma2 * P, XTX)))
    }
    if (verbose && (i %% 100L == 0L || i == samples)) utils::setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)

  list(bet_store   = bet_store,
       lam_store   = lam_store,
       sig_store   = sig_store,
       zeta_store  = zeta_store,
       ed_store    = ed_store,
       accept_rate = acc_post / (samples - burnin),
       sd_prop     = sd_prop,
       samples     = samples,
       burnin      = burnin,
       thin        = thin)
}
