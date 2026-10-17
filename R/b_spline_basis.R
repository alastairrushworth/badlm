# Equally spaced B-spline basis (Eilers and Marx, 1996) evaluated at `x`.
# With `nseg` segments of degree `deg` the basis has `nseg + deg` columns.
b_spline_basis <- function(x, nseg = 10, deg = 3, range.variables = range(x)) {
  x_min <- min(range.variables)
  x_max <- max(range.variables)
  dx    <- (x_max - x_min) / nseg
  knots <- x_min + seq(-deg, nseg + deg) * dx
  splines::splineDesign(knots = knots, x = x, ord = deg + 1, outer.ok = TRUE)
}
