# this one models the prior mean of the last smoothing parameter

gaussian_badlm <- function(X, y, samples, burnin = 0, thin = 1){

  n.save        <- ifelse(thin == 1, (samples - burnin), (samples - burnin) / thin)

  ncolX  <- ncol(X)
  nrowX  <- nrow(X)
  deaths <- y
  XTX    <- crossprod(X)
  D      <- diff(diag(ncolX), diff = 1)
  P      <- crossprod(D)
  Xy     <- t(X) %*% deaths
  yX     <- t(deaths) %*% X

  sig_store  <- vector("numeric", length = n.save)
  max_store  <- vector("numeric", length = n.save)
  eta_store  <- vector("numeric", length = n.save)
  gam_store  <- matrix(0, ncol = ncolX - 1, nrow = n.save)
  bet_store  <- matrix(0, ncol = ncolX, nrow = n.save)
  beta       <- rnorm(ncolX)
  beta[length(beta)] <- 0
  sigma      <- 1
  eta.sq     <- 0.01
  sd.prop    <- 0.1
  z          <- rtrunc(ncolX - 1, mean = 10, sd = 2, a = -10, b = 10, spec = "norm")
  g          <- g_prop <- exp(z)
  Dstar      <- D
  for(k in 1:(ncolX - 1)) Dstar[k, ] <- sqrt(g[k]) * D[k, ]
  P          <- crossprod(Dstar)
  P[nrow(P), nrow(P)] <- P[nrow(P), nrow(P)]
  Q          <- crossprod(diff(diag(ncolX - 1), diff = 1))
  a.P        <- rep(0, ncolX - 1)
  sd_prop    <- rep(1, ncolX - 1)
  tune_per   <- 200
  det.P      <- as.numeric(determinant(P, logarithm = T)$modulus)
  Hat        <- vector()

  ## Start timer
  cat("Collecting", samples, "samples\n", sep = " ")
  progressBar            <- txtProgressBar(style = 3)
  percentage.points      <- round((1:100/100)*samples)
  increment              <- 0
  for(i in 1:samples){
    save.iter <- i > burnin && ((i %% thin == 0) | thin == 0)
    if(save.iter) increment <- increment+1
    # print progress and update the proposal variance
    if(i %% tune_per == 0 ){
      newhat      <- try(sum(diag(X %*% solve(P * sigma  + XTX) %*% t(X))), silent = T)
      if(is.numeric(newhat)) Hat         <- c(Hat, newhat)
      # print(round(z, 0))
      if(i < (samples / 2)){
        prop_change   <- function(p, psd){
          new.sd <- psd
          new.sd <- c(1, 2)[(p > 0.4) + 1]*new.sd
          new.sd <- c(1, 0.5)[(p < 0.3) + 1]*new.sd
          new.sd
        }
      sd_prop        <- prop_change(p = a.P / tune_per, psd = sd_prop)
      a.P            <- a.P*0
      }
    }
    # update beta
    precis         <- (P  + (XTX / sigma))
    beta           <- as.numeric(rmvnorm.canonical(1, b = as.numeric(yX/sigma), Q = precis))

    # update sigma
    a      <- 0.001 + nrowX/2
    b      <- 0.001 + sum((deaths - (X %*% beta))^2)/2
    sigma  <- 1/rgamma(1, shape = a, scale = 1/b)

    # update eta
    a              <- 1 + (ncolX - 1)/2
    b              <- 1 + xQx1(Q, z)/2
    eta.sq         <- 1/rgamma(1, shape = a, scale = 1/b)

    # update P using one-at-a-time
    for(j in 1:length(z)){
      z_prop           <- z
      z_prop[j]        <- rtrunc(1, mean = z_prop[j], sd = sd_prop[j], spec = "norm", b = 40)
      g_prop           <- exp(z_prop)
      Dstar.prop       <- Dstar
      Dstar.prop[j, ]  <- sqrt(g_prop[j]) * D[j, ]
      P_prop                                 <- crossprod(Dstar.prop)
      P_prop[nrow(P_prop), nrow(P_prop)]     <- P_prop[nrow(P_prop), nrow(P_prop)]
      # P_prop              <- P_prop + ridge
      # get the log-likelihood under P_current
      det.P.new        <- as.numeric(determinant(P_prop)$modulus)
      ll_current       <- 0.5*det.P     - (1/2)* xQx1(P, beta)      - (1/(2*eta.sq)) * xQx1(Q, z)
      ll_proposal      <- 0.5*det.P.new - (1/2)* xQx1(P_prop, beta) - (1/(2*eta.sq)) * xQx1(Q, z_prop)
      acceptance       <- exp(ll_proposal - ll_current)
      if(runif(1) <= acceptance){
        z            <- z_prop
        g            <- g_prop
        P            <- P_prop
        det.P        <- det.P.new
        Dstar        <- Dstar.prop
        a.P[j]       <- a.P[j] + 1
      }
    }


    if(save.iter){
      gam_store[increment, ] <- g
      bet_store[increment, ] <- beta
      sig_store[increment]   <- sigma
      eta_store[increment]   <- eta.sq
    }
    if(i %in% percentage.points) setTxtProgressBar(progressBar, i/samples)
  }
  close(progressBar)
  list(sig_store = sig_store, bet_store = bet_store, eta_store = eta_store,
       max_store = max_store, gam_store = gam_store, Hat = Hat)
}
