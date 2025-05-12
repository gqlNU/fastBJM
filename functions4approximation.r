########################################
####  fixed effects
########################################
f_beta_approximated <- function(bx,bx0,params,dat,which_par='') {
  ats <- dat$transitions_to_analyse
  QD <- dat$Q
  pm0 <- params
  nms <- names(bx0)
  pm0[['beta']][nms] <- bx0[nms]
  status <- dat$status
  dev1 <- dev2 <- rep(0,length(bx0))
  names(dev1) <- names(dev2) <- nms
  
  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  hr <- as.numeric(tapply(hq,dat$Q$row_id,sum))
  for (k in ats) {
    tv <- dat$trans_indicators[,k]
    for (xc in dat$XVars) {
      pp <- paste0(xc,'_',k)
      xt <- dat$X[,xc]
      #  first and second derivatives of the full conditionals calculated at bx0
      dvll <- sum(status*tv*xt) - sum(hr*tv*xt)
      dev1[pp] <- dvll - 1/100^2*bx0[pp]
      dvll <- - sum(hr*tv*xt^2)
      dev2[pp] <- dvll - 1/100^2
    }
  }
  out <- approx_coefficients(bx,bx0,dev1,dev2)
  out$x <- c(sapply(1:length(bx0),function(xx) {rnorm(1,out$mean[xx],out$sd[xx])}))
  return(out)
}

########################################
####  random effects in longitudinal
########################################
f_lmm_random_effects_approximated <- function(bx,bx0,params,dat,which_lmm_part) {
  status <- dat$status
  t_ms_t0 <- t_ms_t02 <- t_ms_t0q <- t_ms_t02q <- 1
  if (which_lmm_part=='intercept') pp <- 'u'
  if (which_lmm_part=='slope') {
    pp <- 'd'
    t_ms_t0 <- dat$t_ms_t0
    t_ms_t02 <- t_ms_t0^2
    t_ms_t0q <- dat$Q$t_ms_t0
    t_ms_t02q <- dat$Q$t_ms_t0^2
  }
  pm0 <- params
  pm0[[pp]] <- bx0
  
  ###  from MSM 
  #  hazard at quadrature points
  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  #  longitudinal measurements at event times
  lv <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,1]
  #  longitudinal measurements at quadrature
  lvq <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,1]
  
  a1 <- pm0[['a_1']][dat$jk_index]
  a1q <- pm0[['a_1']][dat$Q$jk_index]
  a2 <- pm0[['a_2']][dat$jk_index]
  a2q <- pm0[['a_2']][dat$Q$jk_index]
  
  #  dev 1
  #    the quadrature part
  ahq <- (a1q*t_ms_t0q + 2*a2q*t_ms_t0q*lvq) * hq
  ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
  #    the event time part
  sat <- status*(a1*t_ms_t0 + 2*a2*t_ms_t0*lv)
  sat_p <- tapply(sat,dat$pid,sum)
  #    putting the two together
  dev1_msm <- sat_p - ahq_p
    
  #  dev 2
  #    the quadrature part
  ahq <- (2*a2q*t_ms_t02q) * hq
  ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
  #    the event time part
  sat <- status*(2*a2*t_ms_t02)
  sat_p <- tapply(sat,dat$pid,sum)
  #    putting the two together
  dev2_msm <- sat_p - ahq_p

  ###  from longitudinal
  dl <- dat$long
  y <- dl$y
  tau_y <- 1/params[['sd_y']]^2
  tt <- rep(1,length(dl$pid))
  if (pp=='d') tt <- dl$t_ms_t0
  #  dev 1
  u <- pm0[['u']][dl$pid]
  d <- pm0[['d']][dl$pid]
  dv <- tau_y*tt*(y - u - d*dl$t_ms_t0)
  dev1_long <- as.numeric(tapply(dv,dl$pid,sum))
  #  dev 2
  dv <- -tau_y*tt^2
  dev2_long <- as.numeric(tapply(dv,dl$pid,sum))
  
  ###  from prior
  tau_RE <- 1/pm0[[paste0('sd_',pp)]]^2
  mean_RE <- pm0[[paste0('m_',pp)]]
  RE <- pm0[[pp]]
  ones <- rep(1,length(RE))
  dev1_prior <- -tau_RE*(RE - mean_RE)
  dev2_prior <- -tau_RE
  
  ###  putting both parts together
  dev1 <- dev1_msm + dev1_long + dev1_prior
  dev2 <- dev2_msm + dev2_long + dev2_prior
  
  out <- approx_coefficients(bx,bx0,dev1,dev2)
  out$x <- c(sapply(1:length(bx0),function(xx) {rnorm(1,out$mean[xx],out$sd[xx])}))
  return(out)
}

########################################
####  association parameters
########################################
f_association_approximated <- function(bx,bx0,params,dat,which_par='') {
  ats <- dat$transitions_to_analyse
  QD <- dat$Q
  pm0 <- params
  nms <- names(bx0)
  pm0[[which_par]][nms] <- bx0[nms]
  status <- dat$status

  dev1 <- dev2 <- rep(0,length(bx0))
  names(dev1) <- names(dev2) <- nms
  
  #  hazards at quadrature points
  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  
  #  linear or quadratic association?
  a_od <- as.numeric(sub('a_','',which_par))
  #  longitudinal values at event times without association parameters
  lh <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,a_od]
  #  longitudinal values at quadrature points without association parameters
  lhq <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,a_od]
  
  d1 <- as.numeric(tapply(hq*lhq,dat$Q$row_id,sum))
  d2 <- as.numeric(tapply(hq*lhq^2,dat$Q$row_id,sum))
  for (k in ats) {
    tv <- dat$trans_indicators[,k]
    dvll <- sum(status*tv*lh) - sum(tv*d1)
    dev1[k] <- dvll - 1/100^2*bx0[k]
    dvll <- - sum(tv*d2)
    dev2[k] <- dvll - 1/100^2
  }
  out <- approx_coefficients(bx,bx0,dev1,dev2)
  out$x <- c(sapply(1:length(bx0),function(xx) {rnorm(1,out$mean[xx],out$sd[xx])}))
  return(out)
}

########################################
####  country-level random effects in MSM
########################################
f_eta_approximated <- function(bx,bx0,params,dat,which_par='') {
  pp <- 'eta'
  ats <- dat$transitions_to_analyse
  QD <- dat$Q
  pm0 <- params
  nms <- names(bx0)
  pm0[[pp]][nms] <- bx0[nms]
  status <- dat$status
  npms <- dat$nctys
  dev1 <- dev2 <- rep(0,length(bx0))
  names(dev1) <- names(dev2) <- nms
  
  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  hr <- as.numeric(tapply(hq,dat$Q$row_id,sum))
  for (k in ats) {
    tv <- dat$trans_indicators[,k]
    for (ix in 1:npms) {
      xt <- dat$cty_dummy[,ix]
      pp <- paste0('cty',ix,'_',k)
      #  first and second derivatives of the full conditionals calculated at bx0
      dvll <- sum(status*tv*xt) - sum(hr*tv*xt)
      dev1[pp] <- dvll - 1/100^2*bx0[pp]
      dvll <- - sum(hr*tv*xt^2)
      dev2[pp] <- dvll - 1/100^2
    }
  }
  out <- approx_coefficients(bx,bx0,dev1,dev2)
  out$x <- c(sapply(1:length(bx0),function(xx) {rnorm(1,out$mean[xx],out$sd[xx])}))
  return(out)
}

########################################
####  country-age random effects in MSM
########################################
f_kappa_approximated <- function(bx,bx0,params,dat,which_par='') {
  pp <- 'kappa'
  ats <- dat$transitions_to_analyse
  QD <- dat$Q
  pm0 <- params
  nms <- names(bx0)
  pm0[[pp]][nms] <- bx0[nms]
  status <- dat$status
  npms <- length(bx0)
  dev1 <- dev2 <- rep(0,length(bx0))
  names(dev1) <- names(dev2) <- nms
  
  #  hazard at quadrature points
  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  #  hazard for each row in the from-to data
  hr <- as.numeric(tapply(hq,dat$Q$row_id,sum))
  cty_age_gps <- c(sapply(1:dat$nGs,function(x){paste0('cty',1:dat$nctys,'_g',x)}))
  for (k in ats) {
    tv <- dat$trans_indicators[,k]
    for (ix in cty_age_gps) {
      pp <- sub('_',paste0('_',k),ix)
      xt <- dat$cty_age_dummy[,ix]
      #  first and second derivatives of the full conditionals calculated at bx0
      dvll <- sum(status*tv*xt) - sum(hr*tv*xt)
      dev1[pp] <- dvll - 1/100^2*bx0[pp]
      dvll <- - sum(hr*tv*xt^2)
      dev2[pp] <- dvll - 1/100^2
    }
  }
  out <- approx_coefficients(bx,bx0,dev1,dev2)
  out$x <- c(sapply(1:length(bx0),function(xx) {rnorm(1,out$mean[xx],out$sd[xx])}))
  return(out)
}


