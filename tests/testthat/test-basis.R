test_that("B-spline basis has k columns, is non-negative and partitions unity", {
  k <- 12
  p <- 30
  B <- badlm:::b_spline_basis(0:p, nseg = k - 3, range.variables = c(0, p))
  expect_equal(dim(B), c(p + 1, k))
  expect_true(all(B >= 0))
  expect_equal(unname(rowSums(B)), rep(1, p + 1))
})
