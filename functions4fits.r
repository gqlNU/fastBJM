
gather_model_spec <- function(model_spec, X_spec) {
  out <- model_spec
  out$byperson <- F
  if (out$msm_by_country=='random') out$include_country_fixed <- F
  out$include_fixed_effects <- F
  if (!is.null(X_spec)) out$include_fixed_effects <- T
  out$joint <- T
  if (out$msm_only) out$joint <- F
  return(out)
}

gather_update_setting <- function(model_spec) {
  out <- list()
  out$update_h0 <- T
  out$update_association <- T
  out$update_fixed_effects <- T
  if (!model_spec$include_fixed_effects) out$update_fixed_effects <- F
  out$update_lmm_random_effects <- T
  out$update_lmm_data_sd <- T
  out$update_cty_random_effects <- T
  out$update_cty_random_effects_SD <- T
  if (model_spec$msm_only) {
    out$update_association <- F
    out$update_lmm_random_effects <- F
    out$update_lmm_data_sd <- F
  }
  if (model_spec$msm_by_country=='fixed' | model_spec$msm_by_country=='') {
    out$update_cty_random_effects <- F
    out$update_cty_random_effects_SD <- F
  }
  out$update_age_cty_interactions <- F
  out$update_age_cty_interactions_SD <- F
  return(out)
}


########################################
####  calculate the log joint posterior
########################################
log_posterior <- function(params, dat) {
  model_spec <- dat$model_spec
  byperson <- model_spec$byperson
  ats <- dat$transitions_to_analyse

  #################################
  #   MSM + LMM
  #################################
  l1 <- ll_msm(params,dat,byperson)
  l2 <- 0
  if (model_spec$joint) l2 <- ll_lmm(params,dat,byperson)
  ll <- l1 + l2

  #################################
  #   prior contributions
  #################################
  prs <- rep(0,14)
  #  age-specific baselines
  prs[1] <- sum(dgamma(params[['l']],shape=0.01,rate=0.01,log=T))
  if (model_spec$include_fixed_effects) prs[2] <- sum(dnorm(params[['beta']],mean=0,sd=100,log=T))
  if (model_spec$msm_by_country=='random') {
    for (k in ats) {
      sigma <- params[['sd_eta']][k]
      #  sum-2-zero country random effects
      prs[3] <- prs[3] + log_density_exchangeable_s2z(params[['eta']][paste0('cty',1:dat$nctys,'_',k)],sigma)
      #  sum-2-zero country random effect SDs
      prs[4] <- prs[4] + dgamma(1/sigma^2,shape=0.01,rate=0.01,log=T)
    }
  }
  if (model_spec$joint) {
    #  means of varying intercepts and slopes
    pp <- 'm_u'
    prs[5] <- dnorm(params[[pp]],0,1000,log=T)
    pp <- 'm_d'
    prs[6] <- dnorm(params[[pp]],0,1000,log=T)
    #  SDs of varying intercepts and slopes
    pp <- 'sd_u'
    prs[7] <- dgamma(1/params[[pp]]^2,shape=0.01,rate=0.01,log=T)
    pp <- 'sd_d'
    prs[8] <- dgamma(1/params[[pp]]^2,shape=0.01,rate=0.01,log=T)
    pp <- 'sd_y'
    prs[9] <- dgamma(1/params[[pp]]^2,shape=0.01,rate=0.01,log=T)
    #  association parameters
    prs[10] <- sum(dnorm(unlist(params[grep('a_',names(params))]),0,100,log=T))
    if (byperson) {
      #  random intercepts
      prs1 <- dnorm(params[['u']],mean=params[['m_u']],sd=params[['sd_u']],log=T)
      #  random slopes
      prs2 <- dnorm(params[['d']],mean=params[['m_d']],sd=params[['sd_d']],log=T)
    } else {
      #  random intercepts
      prs[11] <- sum(dnorm(params[['u']],mean=params[['m_u']],sd=params[['sd_u']],log=T))
      #  random slopes
      prs[12] <- sum(dnorm(params[['d']],mean=params[['m_d']],sd=params[['sd_d']],log=T))
    }
  }
  if (model_spec$msm_country_age) {
    for (k in ats) {
      sigma <- params[['sd_kappa']]
      pms <- c(sapply(paste0('cty',1:dat$nctys,'_',k),function(xx){paste0(xx,'g',1:dat$nGs)}))
      #  sum-2-zero country random effects
      prs[13] <- prs[13] + log_density_exchangeable_s2z(params[['kappa']][pms],sigma)
      #  sum-2-zero country random effect SDs
      prs[14] <- prs[14] + dgamma(1/sigma^2,shape=0.01,rate=0.01,log=T)
    }
  }
  out <- ll + sum(prs)
  if (byperson) out <- out + prs1 + prs2
  return(out)
}
if (F) {
  log_posterior(params,dat)
}

########################################
####  prepare parameters
########################################
get_parameters <- function(dat,model_spec) {
  ats <- dat$transitions_to_analyse
  ps <- c(sapply(ats,function(k){paste0(k,'g',1:dat$nGs)}))
  params <- list(l=rep(0.05,length(ps)))
  names(params[['l']]) <- ps
  if (model_spec$msm_by_country=='random') {
    params[['eta']] <- rep(0,dat$nctys*dat$NTS)
    names(params[['eta']]) <- sapply(ats,function(x){paste0('cty',1:dat$nctys,'_',x)})
    params[['sd_eta']] <- rep(0.1,dat$NTS)
    names(params[['sd_eta']]) <- ats
  }
  if (model_spec$include_fixed_effects) {
    params[['beta']] <- rep(0,dat$nXs*dat$NTS)
    names(params[['beta']]) <- sapply(ats,function(k){paste0(dat$XVars,'_',k)})
  }
  if (model_spec$msm_country_age) {
    params[['kappa']] <- rep(0,dat$nctys*dat$NTS*dat$nGs)
    names(params[['kappa']]) <- sapply(sapply(ats,function(x){paste0('cty',1:dat$nctys,'_',x)}),function(xx){paste0(xx,'g',1:dat$nGs)})
    params[['sd_kappa']] <- 0.1
  }
  if (model_spec$joint) {
    params[['a_1']] <- rep(0,dat$NTS)
    names(params[['a_1']]) <- ats
    params[['a_2']] <- params[['a_1']]
    params[['u']] <- params[['d']] <- rep(0,dat$npds)
    params[['sd_y']] <- 1
    params[['sd_u']] <- 1
    params[['sd_d']] <- .1
    params[['m_u']] <- 0
    params[['m_d']] <- 0
  }
  return(params)
}


#######################################################
###  function to calculate the fixed effects
###  on the log scale
#######################################################
compute_fixed_effects <- function(params, dat, at_quadrature) {
  X <- dat$X
  ats <- dat$transitions_to_analyse
  trans_indicators <- dat$trans_indicators
  bX <- NULL
  for (k in ats) {
    #  extract beta
    beta <- matrix(params[['beta']][paste0(dat$XVars,'_',k)],ncol=1)
    #  calculate X * beta for that transition
    bX <- cbind(bX,X%*%beta*trans_indicators[,k])
  }
  out <- apply(bX,1,sum)
  if (at_quadrature & dat$nQs>1) out <- rep(out,each=dat$nQs)
  return(out)
}
if (F) {
  summary(compute_fixed_effects(params,all_fitdata,F))
}

#######################################################
###  function to calculate longitudinal effects
#######################################################
compute_longitudinal_effects <- function(params, dat, at_quadrature, include_association_par=T) {
  #  get the relevant data
  dd <- dat
  if (at_quadrature) dd <- dat$Q
  pid <- dd$pid
  jk_index <- dd$jk_index
  t_ms_t0 <- dd$t_ms_t0
  
  #  extract and format parameters
  #  a1 and a2 as in a1*BMI and a2*BMI^2
  A <- cbind(params[['a_1']][jk_index],
             params[['a_2']][jk_index])
  #  compute (standardised) longitudinal values
  u <- params[['u']][pid]
  d <- params[['d']][pid]
  sv <- u + d*t_ms_t0
  sv2 <- sv*sv
  #  BMI , BMI^2
  out <- cbind(sv,sv2)
  if (include_association_par) out <- rowSums(A*out)
  return(out)
}

#######################################################
###  function to calculate the hazard ratios
#######################################################
compute_RR <- function(params, dat, at_quadrature) {
  model_spec <- dat$model_spec
  dd <- dat
  if (at_quadrature) dd <- dat$Q
  bX <- aL <- eta <- kappa <- 0
  #  fixed effects
  if (model_spec$include_fixed_effects) bX <- compute_fixed_effects(params, dat, at_quadrature)
  #  longitudinal
  if (model_spec$joint) aL <- compute_longitudinal_effects(params, dat, at_quadrature, include_association_par=T)
  #  country-level random effects
  if (model_spec$msm_by_country=='random') eta <- params[['eta']][dd$jkc_index]
  #  age-country random effects
  if (model_spec$msm_country_age) kappa <- params[['kappa']][dd$jkgc_index]
  #  putting all together
  logRR <- bX + aL + eta + kappa
  #  return hazard ratios
  out <- exp(logRR)
  return(out)
}

#######################################################
###  function to calculate the transition hazards
#######################################################
compute_haz <- function(params, dat, at_quadrature, include_quadweights=T, without='') {
  h0 <- compute_baseline(params, dat, at_quadrature)
  RR <- compute_RR(params, dat, at_quadrature)
  out <- h0*RR
  if (without=='h0') out <- RR
  if (at_quadrature & include_quadweights) out <- out * dat$Q$w  #  multiplying quadrature weights
  return(out)
}

#######################################################
###  function to calculate the multistate likelihood
#######################################################
ll_msm <- function(params, dat, byperson=FALSE) {
  status <- dat[['status']]
  ll1 <- compute_haz(params,dat,at_quadrature = F)^status
  ll2 <- compute_haz(params,dat,at_quadrature = T)  # this includes quadrature weights
  if (byperson) {
    tll1 <- tapply(log(ll1),dat$pid,sum)
    tll2 <- tapply(ll2,dat$Q$pid,sum)
    out <- tll1 - tll2
  } else {
    out <- sum(log(ll1)) - sum(ll2)
  }
  return(out)
}
if (F) {
  dat <- all_fitdata
  ll_msm(params,dat)
}


#######################################################
###  function to calculate the longitudinal likelihood
#######################################################
ll_lmm <- function(params, dat, byperson=F) {
  dd <- dat$long
  y <- dd$y
  tt <- dd$t_ms_t0
  pid <- dd$pid
  u <- params[['u']][pid] # intercepts
  d <- params[['d']][pid] # slopes
  m <- u + d*tt
  sy <- params[['sd_y']]
  ll <- dnorm(y,mean=m,sd=sy,log=T)
  if (byperson) {
    out <- tapply(ll,pid,sum)
  } else {
    out <- sum(ll)
  }
  return(out)
}

#######################################################
###  function to update age-specific baselines
#######################################################
update_baseline <- function(params, dat, prior=c(0.01,0.01), reference_ids=NULL) {
  nms <- names(params[['l']])
  shape <- rep(0,length(nms))
  names(shape) <- nms
  if (!is.null(reference_ids)) {
    ii <- reference_ids
    ts <- tapply(dat$status[ii],dat$jkg_index[ii],sum)
  } else {
    ts <- dat$ncases_by_jkg
  }
  shape[nms] <- ts[nms]
  wh <- compute_haz(params, dat, at_quadrature = T, without='h0')
  whc <- tapply(wh,dat$Q$jkg_index,sum)
  rate <- rep(0,length(nms))
  names(rate) <- nms
  rate[names(whc)] <- whc
  out <- sapply(1:length(shape),function(x){
       rgamma(1,shape=shape[x]+prior[1],rate=rate[x]+prior[2])})
  if (any(out<=0)) out[which(out<=0)] <- 1e-30
  names(out) <- nms
  return(out)
}

#######################################################
###  function to find the mode iteratively
#######################################################
auto_mode_finder <- function(which_par, start_x0, params, dat) {
  # updated on 31Jan2025: fixed the problem that v2_31Jan2025 returns 
  #                       a value found at the penultimate step, not the 
  #                       final optimal one
  gmrf_settings <- dat$gmrf_settings
  nm <- names(start_x0)
  x0 <- start_x0
  iter <- 0
  xr <- gmrf_settings$search_range
  sp <- gmrf_settings$search_points
  xv <- seq(xr[1],xr[2],length.out=sp)

  if (length(grep('a_',which_par))>0) fa <- f_association_approximated
  if (which_par=='beta') fa <- f_beta_approximated
  if (which_par=='intercept' | which_par=='slope') fa <- f_lmm_random_effects_approximated
  if (which_par=='eta') fa <- f_eta_approximated
  if (which_par=='kappa') fa <- f_kappa_approximated

  d <- fa(xv,x0,params,dat,which_par)$d
  new_x0 <- apply(d,2,function(xx){xv[which.max(xx)]})
  names(new_x0) <- nm
  iter <- iter + 1
  err <- abs((x0-new_x0)/new_x0)*100
  while (any(err>gmrf_settings$eps) & iter < gmrf_settings$max_iters) {
    x0 <- new_x0
    d <- fa(xv,x0,params,dat,which_par)$d
    new_x0 <- apply(d,2,function(xx){xv[which.max(xx)]})
    names(new_x0) <- nm
    err <- abs((x0-new_x0)/new_x0)*100
    iter <- iter + 1
  }
  out <- list(x0=new_x0,iter=iter,abs_err100=err)
  return(out)
}

#######################################################
###  function to sample parameter values based on GMRF
#######################################################
gmrf_sampling <- function(which_par, params, dat) {
  if (which_par=='intercept') pp <- 'u'
  if (which_par=='slope') pp <- 'd'
  if (length(grep('a_',which_par))>0) pp <- which_par
  if (which_par=='beta') pp <- which_par
  start_x0 <- params[[pp]]
  
  opt_x <- auto_mode_finder(which_par,start_x0,params,dat)
  if (which_par=='intercept' | which_par=='slope') fa <- f_lmm_random_effects_approximated
  if (length(grep('a_',which_par))>0) fa <- f_association_approximated
  if (which_par=='beta') fa <- f_beta_approximated
  apx <- fa(0,opt_x$x0,params,dat,which_par) # first argument does not matter
  names(apx$x) <- names(start_x0)
  out <- list(mean=apx$mean,sd=apx$sd,x=apx$x)
  return(out)
}

#######################################################
###  function to calculate the coefficients in the
###  approximation aa+bb*x-1/2cc*x^2
###  with aa, bb and cc calculated at x0
###
###  aa <- f0 - x0*dev1 + x0^2/2*dev2
#######################################################
approx_coefficients <- function(x,x0,dev1,dev2) {
  bb <- dev1 - x0*dev2
  cc <- -dev2
  mm <- bb/cc  # mean
  prec <- cc   # precision
  #  approximated density at x using Taylor expansion at x0
  d <- matrix(0,nrow=length(x),ncol=length(x0))
  #for (ix in 1:length(x0)) d[,ix] <- f0[ix] + dev1[ix]*(x-x0[ix]) + dev2[ix]/2*(x-x0[ix])^2
  for (ix in 1:length(x0)) d[,ix] <- dev1[ix]*(x-x0[ix]) + dev2[ix]/2*(x-x0[ix])^2
  #  putting everything together for output
  out <- list(b=bb,c=cc,d=d)
  out$mean <- mm
  out$sd <- sqrt(1/prec)
  return(out)
}

update_random_effect_mean <- function(params,dat,which_lmm_part) {
  #
  #  prior is in the form of N(prior[1]=0,prior[2]=precision)
  #
  # 	which_random_effect <- 'u'  #  intercept
  # 	which_random_effect <- 'd'  #  slope
  #
  if (which_lmm_part=='intercept') pp <- 'u'
  if (which_lmm_part=='slope') pp <- 'd'
  
  tau_RE <- 1/params[[paste0('sd_',pp)]]^2
  tau <- 1/100^2
  
  N <- dat$npds
  RE <- params[[pp]]
  
  m <- (tau_RE*N)/(tau_RE*N + tau) * mean(RE)
  p <- tau_RE*N + tau
  out <- rnorm(1,mean=m,sd=1/sqrt(p))
  return(out)
}
if (F) {
  update_random_effect_mean(inits1,final,'intercept')
}

update_random_effect_sd <- function(params,dat,which_lmm_part) {
  #
  # which_random_effect <- 'u'  #  intercept
  # which_random_effect <- 'b'  #  slope
  #
  if (which_lmm_part=='intercept') pp <- 'u'
  if (which_lmm_part=='slope') pp <- 'd'
  
  m_RE <- params[[paste0('m_',pp)]]
  RE <- params[[pp]]
  N <- dat$npds
  
  shape <- N/2
  prior <- c(0.01,0.01)
  rate <- 0.5*sum((RE-m_RE)^2)
  prec <- rgamma(1,shape=shape+prior[1],rate=rate+prior[2])
  out <- sqrt(1/prec)  #  return SD
  return(out)
}
if (F) {
  update_random_effect_sd(inits1,final, which_lmm_part ='intercept')
}

update_data_sd <- function(params,dat) {
  prior <- c(0.01,0.01)
  dl <- dat$long
  y <- dl$y      #  longitudinal observations
  tt <- dl$t_ms_t0
  nts <- dl$nts  #  number of observations per person
  pid <- dl$pid  #  which y belongs to which person
  u <- params[['u']][pid] # intercept
  d <- params[['d']][pid] # slope
  shape <- sum(nts)/2
  rate <- 0.5*sum((y-u-d*tt)^2)
  prec <- rgamma(1,shape=shape+prior[1],rate=rate+prior[2])
  #  returnn sd
  out <- sqrt(1/prec)
  return(out)
}
if (F) {
  update_data_sd(inits1,final)
}


sample_country_random_effects <- function(params, dat, eval=NULL) {
  #  based on the method on p. 171 in Rue and Held 2005
  all_can_app <- canonical_approximation('eta',params,dat)
  ats <- dat$transitions_to_analyse
  out <- as.list(1:dat$NTS)
  names(out) <- ats
  for (k in ats) {
    can_app <- all_can_app[[k]]  # for the k(th) transition
    A <- can_app$A
    #  change the canonical parameterisation to the usual normal distribution
    norm_app <- canonical_to_normal(can_app$Q_mu,can_app$Q)
    if (is.null(eval)) {
      ###################################################
      #   sample then evaluate p(proposed | current)
      ###################################################
      #  sample from the normal with a sum to zero constrant (i.e. Algorithm 2.6)
      #    output has both the constrained (xstar) and unconstrained (x) sample 
      sm <- algorithm_2_6(mu=norm_app$mu,Cov=norm_app$Cov)
      x <- matrix(sm$x,nrow=dat$nctys)
      Ax <- matrix(sm$xstar,nrow=dat$nctys)
    } else {
      ###################################################
      #   evaluate p(current | proposed)
      ###################################################
      x <- eval$x[paste0('cty',1:dat$nctys,'_',k)]
      Ax <- eval$Ax[paste0('cty',1:dat$nctys,'_',k)]
    }
    #  log density of the unconstrained x  (p(x))
    ld1 <- mvtnorm::dmvnorm(c(x),norm_app$mu,sigma=norm_app$Cov,log=T)
    #  log density of the constrained sample given the unconstrained one log(p(Ax|x))
    ld2 <- -1/2*log(det(A%*%t(A)))
    #  log density of the constraint sample Ax
    m <- A%*%matrix(norm_app$mu,nrow=dat$nctys)
    Cv <- A%*%t(MASS::ginv(norm_app$Cov))%*%t(A)
    ld3 <- dnorm(A%*%Ax,mean=m,sd=sqrt(Cv),log=T)
    ld <- ld1 + ld2 - ld3
    tout <- list(x=x,Ax=Ax,log_den=ld)
    if (!is.null(eval)) tout <- list(log_den=ld)
    out[[k]] <- tout
  }
  if (is.null(eval)) {
    nms <- tx <- tAx <- tld <- NULL
    for (k in ats) {
      nms <- c(nms,paste0('cty',1:dat$nctys,'_',k))
      tx <- c(tx,out[[k]]$x)
      tAx <- c(tAx,out[[k]]$Ax)
      tld <- c(tld,out[[k]]$log_den)
    }
    names(tld) <- ats
    names(tx) <- names(tAx) <- nms
    out <- list(x=tx,Ax=tAx,log_den=tld)
  } else {
    tld <- NULL
    for (k in ats) tld <- c(tld,out[[k]]$log_den)
    names(tld) <- ats
    out <- list(log_den=tld)
  }
  return(out)
}
if (F) {
  out <- sample_country_random_effects(params=params,dat=all_fitdata)
}

sample_age_country_interactions <- function(params, dat, eval=NULL) {
  #  based on the method on p. 171 in Rue and Held 2005
  all_can_app <- canonical_approximation('kappa',params,dat)
  ats <- dat$transitions_to_analyse
  out <- as.list(1:dat$NTS)
  names(out) <- ats
  for (k in ats) {
    can_app <- all_can_app[[k]]  # for the k(th) transition
    A <- can_app$A
    #  change the canonical parameterisation to the usual normal distribution
    norm_app <- canonical_to_normal(can_app$Q_mu,can_app$Q)
    if (is.null(eval)) {
      ###################################################
      #   sample then evaluate p(proposed | current)
      ###################################################
      #  sample from the normal with a sum to zero constrant (i.e. Algorithm 2.6)
      #    output has both the constrained (xstar) and unconstrained (x) sample 
      sm <- algorithm_2_6(mu=norm_app$mu,Cov=norm_app$Cov)
      x <- matrix(sm$x,nrow=dat$nctys*dat$nGs)
      Ax <- matrix(sm$xstar,nrow=dat$nctys*dat$nGs)
    } else {
      ###################################################
      #   evaluate p(current | proposed)
      ###################################################
      pms <- c(sapply(paste0('cty',1:dat$nctys,'_',k),function(xx){paste0(xx,'g',1:dat$nGs)}))
      x <- eval$x[pms]
      Ax <- eval$Ax[pms]
    }
    #  log density of the unconstrained x  (p(x))
    ld1 <- mvtnorm::dmvnorm(c(x),norm_app$mu,sigma=norm_app$Cov,log=T)
    #  log density of the constrained sample given the unconstrained one log(p(Ax|x))
    ld2 <- -1/2*log(det(A%*%t(A)))
    #  log density of the constraint sample Ax
    m <- A%*%matrix(norm_app$mu,nrow=dat$nctys*dat$nGs)
    Cv <- A%*%t(MASS::ginv(norm_app$Cov))%*%t(A)
    ld3 <- dnorm(A%*%Ax,mean=m,sd=sqrt(Cv),log=T)
    ld <- ld1 + ld2 - ld3
    tout <- list(x=x,Ax=Ax,log_den=ld)
    if (!is.null(eval)) tout <- list(log_den=ld)
    out[[k]] <- tout
  }
  if (is.null(eval)) {
    nms <- tx <- tAx <- tld <- NULL
    for (k in ats) {
      pms <- c(sapply(paste0('cty',1:dat$nctys,'_',k),function(xx){paste0(xx,'g',1:dat$nGs)}))
      nms <- c(nms,pms)
      tx <- c(tx,out[[k]]$x)
      tAx <- c(tAx,out[[k]]$Ax)
      tld <- c(tld,out[[k]]$log_den)
    }
    names(tld) <- ats
    names(tx) <- names(tAx) <- nms
    out <- list(x=tx,Ax=tAx,log_den=tld)
  } else {
    tld <- NULL
    for (k in ats) tld <- c(tld,out[[k]]$log_den)
    names(tld) <- ats
    out <- list(log_den=tld)
  }
  return(out)
}
if (F) {
  out <- sample_age_country_interactions(params=params,dat=all_fitdata)
}

