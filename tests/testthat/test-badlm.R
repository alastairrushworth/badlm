sim_dlm <- function(n, truth, sd = 0.01) {
  p <- length(truth) - 1
  x <- as.numeric(arima.sim(list(ar = 0.5), n = n))
  y <- as.numeric(lag_matrix(x, p) %*% truth + rnorm(n - p, sd = sd))
  list(x = x, y = y)
}

test_that("badlm recovers a known lag curve", {
  skip_on_cran()
  set.seed(3)
  p     <- 20
  truth <- dnorm(0:p, mean = 6, sd = 3) / 10
  d     <- sim_dlm(400, truth)
  fit   <- badlm(d$x, d$y, k = 10, nlag = p, samples = 2000, verbose = FALSE)
  expect_equal(dim(fit$coefs), c(p + 1, 4))
  expect_equal(dim(fit$lag_samples), c(1000, p + 1))
  expect_length(fit$ed, 1000)
  expect_true(all(fit$coefs$lower <= fit$coefs$beta & fit$coefs$beta <= fit$coefs$upper))
  rmse <- sqrt(mean((fit$coefs$beta - truth)^2))
  expect_lt(rmse, 0.15 * max(truth))
  coverage <- mean(truth >= fit$coefs$lower & truth <= fit$coefs$upper)
  expect_gt(coverage, 0.7)
})

test_that("badlm accepts a pre-truncated response and validates its inputs", {
  set.seed(4)
  x <- rnorm(50)
  y <- rnorm(50)
  f_full <- badlm(x, y, k = 5, nlag = 5, samples = 60, burnin = 20, verbose = FALSE)
  f_trim <- badlm(x, y[-(1:5)], k = 5, nlag = 5, samples = 60, burnin = 20, verbose = FALSE)
  expect_equal(nrow(f_full$coefs), 6)
  expect_equal(nrow(f_trim$coefs), 6)
  expect_equal(nrow(f_full$lag_samples), 40)
  expect_error(badlm(x, y[1:30], k = 5, nlag = 5, samples = 60, verbose = FALSE))
  expect_error(badlm(x, y, k = 3, nlag = 5, samples = 60, verbose = FALSE))
  expect_error(badlm(c(x, NA), c(y, 1), k = 5, nlag = 5, samples = 60, verbose = FALSE))
})

test_that("plot_lagcurve returns a ggplot", {
  set.seed(5)
  x <- rnorm(50)
  y <- rnorm(50)
  fit <- badlm(x, y, k = 5, nlag = 5, samples = 60, burnin = 20, verbose = FALSE)
  expect_s3_class(plot_lagcurve(fit), "ggplot")
})
