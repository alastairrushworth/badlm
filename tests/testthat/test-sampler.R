test_that("the adaptive penalty is singular with the closed-form generalised determinant", {
  # The sampler uses log pdet(P) = log(K) + sum(tau) for P = D' diag(exp(tau)) D.
  set.seed(1)
  K <- 12
  D <- diff(diag(K))
  for (s in 1:5) {
    tau <- rnorm(K - 1, mean = 5, sd = 2)
    P   <- crossprod(sqrt(exp(tau)) * D)
    ev  <- sort(eigen(P, symmetric = TRUE, only.values = TRUE)$values)
    expect_lt(abs(ev[1]), 1e-8 * max(ev))
    expect_equal(sum(log(ev[-1])), log(K) + sum(tau))
  }
})

test_that("gaussian_badlm returns draws of the right shape and honours burn-in and thinning", {
  set.seed(2)
  n <- 60
  K <- 6
  X <- matrix(rnorm(n * K), n, K)
  y <- as.numeric(X %*% rep(0.5, K) + rnorm(n))
  fit <- badlm:::gaussian_badlm(X, y, samples = 300, burnin = 100, thin = 2, verbose = FALSE)
  expect_equal(dim(fit$bet_store), c(100, K))
  expect_equal(dim(fit$lam_store), c(100, K - 1))
  expect_length(fit$sig_store, 100)
  expect_length(fit$zeta_store, 100)
  expect_length(fit$ed_store, 100)
  expect_true(all(is.finite(fit$bet_store)))
  expect_true(all(fit$sig_store > 0))
  expect_true(all(fit$ed_store > 0 & fit$ed_store <= K))
  expect_true(all(fit$accept_rate >= 0 & fit$accept_rate <= 1))
  expect_error(badlm:::gaussian_badlm(X, y, samples = 100, burnin = 100))
})
