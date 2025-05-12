
check_all_elements_the_same <- function(x,eps) {
  out <- T
  d <- abs((x-x[1])/x[1])
  if (any(d>eps)) out <- F
  return(out)
}
IGMRF_Q <- function(n,s) {
  #  s = standard deviation
  vsmall <- 1e-10
  Q <- array(vsmall,c(n,n))
  tau <- 1/s^2
  diag(Q) <- tau
  dg <- eigen(Q)
  val <- dg$values
  vec <- dg$vectors
  e1 <- vec[,1]
  if (!check_all_elements_the_same(e1,0.05)) {
    print(e1)
    print(dg)
    stop('something strange in the eigenvectors')
  }
  b <- c(0,val[2:n])
  tQ <- vec%*%diag(b)%*%t(vec)
  A <- t(e1)
  out <- list(tQ=tQ,A=A,val=val)
  return(out)
}

log_density_exchangeable_s2z <- function(x,sd_x) {
  #  this is the density of a set of exchangeable normal random effects
  #  with the sum to zero constraint (p. 90 in Rue and Held 2005)
  n <- length(x)
  iQ <- IGMRF_Q(n,sd_x)
  tQ <- iQ$tQ
  val <- iQ$val
  out <- c(-(n-1)/2*log(2*pi) + 0.5*sum(log(val[2:n])) - 0.5*matrix(x,nrow=1)%*%tQ%*%matrix(x,nrow=n))
  return(out)
}

canonical_approximation <- function(which_par, params, dat) {
  #  based on the method on p. 171 in Rue and Held 2005
  out_b <- out_c <- NULL
  pp <- which_par #'eta'
  v0 <- params[[pp]]
  opt_v0 <- iters <- NULL
  
  x0c <- auto_mode_finder(pp, v0, params, dat)
  opt_v0 <- x0c$x0
  if (pp=='eta') app <- f_eta_approximated(0,x0c$x0,params,dat) # no need to worry about the first argument
  if (pp=='kappa') app <- f_kappa_approximated(0,x0c$x0,params,dat) # no need to worry about the first argument
  out_b <- app$b
  out_c <- app$c
  
  ats <- dat$transitions_to_analyse
  out <- as.list(1:dat$NTS)
  names(out) <- ats
  for (k in ats) {
    if (pp=='eta') {
      sigma <- params[['sd_eta']][k]
      pms <- paste0('cty',1:dat$nctys,'_',k)
    }
    if (pp=='kappa') {
      sigma <- params[['sd_kappa']]
      pms <- c(sapply(paste0('cty',1:dat$nctys,'_',k),function(xx){paste0(xx,'g',1:dat$nGs)}))
    }
    bk <- out_b[pms]
    ck <- out_c[pms]
    B <- matrix(bk,nrow=length(bk))
    C <- ck
    gQ <- IGMRF_Q(length(bk),sigma)
    tQ <- gQ$tQ
    A <- gQ$A  #  the sum2zero constraint applied to x
    mu <- matrix(0,nrow=length(bk))
    # under canonical parameterisation NC(Qmu,Q)
    Q_mu <- tQ%*%mu + B
    Q <- tQ + diag(C)
    out[[k]] <- list(Q_mu=Q_mu,Q=Q,A=A)
  }
  return(out)
}

canonical_to_normal <- function(Q_mu,Q) {
  # NC(Qmu,Q) = N(mu,Qinv)
  Qinv <- MASS::ginv(Q)
  mu <- Qinv%*%Q_mu
  out <- list(mu=mu,Cov=Qinv)
  return(out)
}

algorithm_2_6 <- function(mu,Cov) {
  #  as on p. 38 in Rue and Held 2005 a the sum to zero constraint
  n <- nrow(Cov)
  Q <- MASS::ginv(Cov)
  L <- chol(Q)
  z <- rnorm(n)
  v <- MASS::ginv(t(L))%*%matrix(z,nrow=n)[,1]
  x <- mu + v
  A <- matrix(rep(1,n),nrow=1)
  V <- Cov%*%t(A)
  W <- A%*%V
  U <- MASS::ginv(W)%*%t(V)
  cc <- A%*%matrix(x,nrow=n) - 0
  xstar <- c(x - t(U)%*%matrix(cc,nrow=1))
  if (abs(sum(xstar))>0.001) stop('values are not summed to 0')
  out <- list(x=c(x),xstar=xstar)
  return(out)
}
if (F) {
  n <- 6
  algorithm_2_6(mu=rep(0,n),Cov=diag(rep(4,n)))
}
if (F) {
  params <- inits1
  dat <- final
  which_transition <- '12'
  can_app <- canonical_approximation(params,dat,which_transition)
  nm <- canonical_to_normal(can_app$Q_mu,can_app$Q)
  algorithm_2_6(mu=nm$mu,Cov=nm$Cov)
}
