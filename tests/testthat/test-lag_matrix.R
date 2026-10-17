# Reference implementation from badlm 0.0.0.9000
lag_matrix_ref <- function(rain, p, start.at.zero = TRUE) {
  slider <- function(n) rev(rain[n:(n + p - !start.at.zero)])
  nums   <- as.list(1:(length(rain) - p))
  matrix(unlist(lapply(nums, slider)), nrow = length(rain) - p,
         ncol = p + start.at.zero, byrow = TRUE)
}

test_that("lag_matrix matches the original implementation", {
  set.seed(1)
  x <- rnorm(40)
  for (p in c(1, 3, 10)) {
    expect_equal(lag_matrix(x, p), lag_matrix_ref(x, p))
    expect_equal(lag_matrix(x, p, start.at.zero = FALSE),
                 lag_matrix_ref(x, p, start.at.zero = FALSE))
  }
})

test_that("lag_matrix has the documented layout and validates p", {
  m <- lag_matrix(1:6, p = 2)
  expect_equal(dim(m), c(4, 3))
  expect_equal(m[1, ], c(3, 2, 1))
  expect_equal(m[4, ], c(6, 5, 4))
  expect_equal(dim(lag_matrix(1:6, p = 2, start.at.zero = FALSE)), c(4, 2))
  expect_error(lag_matrix(1:6, p = 6))
  expect_error(lag_matrix(1:6, p = 0))
})
