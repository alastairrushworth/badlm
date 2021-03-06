#' @export
#' @importFrom Rcpp compileAttributes
badlm <- function(x, y, k, nlag, samples = 10000, thin = 1){
  # truncate the last nlag values of y
  if(length(x) == length(y)) y <- y[(nlag + 1):length(y)]
  
  # construct the lag matrix
  lag_mat       <- lag_matrix(x, p = nlag)

  # transform the lag_mat through a B-spline basis
  bbase         <- b_spline_basis(0:nlag, nseg = k - 3, range.variables = range(0:nlag))
  X_mat         <- rbind(as.matrix(lag_mat %*% bbase), bbase[nlag + 1,])
  
  # fit the distributed lag model
  model_out    <- gaussian_badlm(X = X_mat, y = c(y, 0), samples = samples, thin = thin)
  
  # return the output
  fcoefs <- as.numeric(bbase %*% apply(model_out$bet_store, 2, quantile, probs = 0.5))
  fupper <- as.numeric(bbase %*% apply(model_out$bet_store, 2, quantile, probs = 0.975))
  flower <- as.numeric(bbase %*% apply(model_out$bet_store, 2, quantile, probs = 0.025))
  ftable <- tibble(lag = 0:nlag, beta = fcoefs, upper = fupper, lower = flower)
  
  return(list(coefs = ftable, model_out))
}
