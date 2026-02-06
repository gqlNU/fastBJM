########################################
####  updating parameters via 
####  Metropolis-within-Gibbs
########################################
#' @export
mcmc_update <- function(data, inits, niters, model_spec, update_setting) {
  data$model_spec$byperson <- F
  current_pars <- inits
  sims.list <- setup_storage(current_pars)
  current_jpd <- log_posterior(current_pars, data)

  ##   initial on screen display
  t1 <- Sys.time()
  msg <- paste0(' ** fitting started at ',t1)
  cat(msg)

  ##   initiate progress bar
  pb_freq <- ifelse(niters<10,1,10)  ##   how often to update pb
  pb <- txtProgressBar(min=0, max=niters, initial=0, style=3)

  for (iter in 1:niters) {
      ##-##################################################
      ##  age-specific baselines
      ##-##################################################
      if (update_setting$update_h0) {
          pp <- 'l'
          x <- update_baseline(current_pars, data, prior=c(0.01,0.01))
          nms <- names(x)
          current_pars[[pp]][nms] <- x[nms]
          sims.list[[pp]] <- rbind(sims.list[[pp]],current_pars[[pp]])
          current_jpd <- log_posterior(current_pars, data)
      }
      ##-#########################################################
      ##   update intercepts and slopes (p. 169 in Rue and Held)
      ##-#########################################################
      if (update_setting$update_lmm_random_effects) {
          u0 <- cbind(current_pars[['u']],current_pars[['d']])
          which_lmm_part <- 'lmm_corr_random_effects'
          ##   update the correlted random effects
          propose_pars <- current_pars
          current_jpd_byperson <- current_jpd
          if (model_spec$update_RE_byperson) {
              ##   update by person
              data$model_spec$byperson <- T
              current_jpd_byperson <- log_posterior(current_pars, data)
          }
          upt <- gmrf_sampling(which_lmm_part,current_pars,data)
          ustar <- upt$x
          propose_pars[['u']] <- ustar[,1] ##   intercepts
          propose_pars[['d']] <- ustar[,2] ##   slopes
          names(propose_pars[['u']]) <- names(propose_pars[['d']]) <- data$subject
          propose_jpd_byperson <- log_posterior(propose_pars,data)
          ##  q(ustar|u0)
          mp <- check_positive_definiteness(upt$V)
          ind_lq_bottom <- sapply(1:data$npds,
                                     function(xx){mvtnorm::dmvnorm(ustar[xx,],
                                                          mean=upt$mean[xx,],
                                                          sigma=round(upt$V[xx,,],digits=10),
                                                          log=T)})
          lq_bottom <- sum(ind_lq_bottom)
          ##  for q(u0|ustar)
          upt_propose <- gmrf_sampling(which_lmm_part,propose_pars,data)
          tmp <- check_positive_definiteness(upt_propose$V)
          ind_lq_top <- sapply(1:data$npds,
                                   function(xx){mvtnorm::dmvnorm(u0[xx,],
                                                        mean=upt_propose$mean[xx,],
                                                        sigma=round(upt_propose$V[xx,,],digits=10),
                                                        log=T)})
          lq_top <- sum(ind_lq_top)
          d <- propose_jpd_byperson - current_jpd_byperson + ind_lq_top - ind_lq_bottom
          if (length(which(d< -99999))>0) d[which(d< -99999)] <- -99999
          if (!model_spec$update_RE_byperson) d <- sum(d)
          lr <- log(runif(length(d)))
          if (model_spec$update_RE_byperson) {
              ##   update each person independently
              icc <- which(d>=lr)
              for (pp in c('u','d')) {
                  current_pars[[pp]][icc] <- propose_pars[[pp]][icc]
              }
              ##   calculate the joint posterior based on the updated parameters
              data$model_spec$byperson <- F
              current_jpd <- log_posterior(current_pars,data)
          } else {
              ##   update all people in one go
              if (d>=lr) {
                  for (pp in c('u','d')) current_pars[[pp]] <- propose_pars[[pp]]
                  current_jpd <- propose_jpd_byperson
              }
          }
          ##   store iteration
          for (pp in c('u','d')) sims.list[[pp]] <- rbind(sims.list[[pp]],current_pars[[pp]])
      }
      ##-##########################
      ##   update fixed effects
      ##-##########################
      if (update_setting$update_fixed_effects) {
          pp <- 'beta'
          propose_pars <- current_pars
          bstar <- gmrf_sampling(pp, current_pars, data)
          propose_pars[[pp]][update_setting$beta_to_update] <- bstar$x[update_setting$beta_to_update]
          ##   q(ustar|u0)
          lq_bottom <- sum(dnorm(bstar$x[update_setting$beta_to_update],
                                 mean=bstar$mean[update_setting$beta_to_update],
                                 sd=bstar$sd[update_setting$beta_to_update],log=T))
          bc <- gmrf_sampling(pp, propose_pars, data)
          ##   q(u0|ustar)
          lq_top <- sum(dnorm(current_pars[[pp]][update_setting$beta_to_update],
                              mean=bc$mean[update_setting$beta_to_update],
                              sd=bc$sd[update_setting$beta_to_update],log=T))
          mh_factor <- lq_top - lq_bottom
          propose_jpd <- log_posterior(propose_pars, data)
          d <- propose_jpd - current_jpd + mh_factor
          if (d>=0) {
              current_pars <- propose_pars
              current_jpd <- propose_jpd
          } else {
              if (d>log(runif(1))) {
                  current_pars <- propose_pars
                  current_jpd <- propose_jpd
              }
          }
          sims.list[[pp]] <- rbind(sims.list[[pp]],current_pars[[pp]])
      }
      ##-################################
      ##   update association parameters
      ##-################################
      if (update_setting$update_association) {
          for (pp in model_spec$aps) {
              propose_pars <- current_pars
              bstar <- gmrf_sampling(pp, current_pars, data)
              propose_pars[[pp]] <- bstar$x
              ##   q(ustar|u0)
              lq_bottom <- dnorm(bstar$x,mean=bstar$mean,sd=bstar$sd,log=T)
              bc <- gmrf_sampling(pp, propose_pars, data)
              ##   q(u0|ustar)
              lq_top <- dnorm(current_pars[[pp]],mean=bc$mean,sd=bc$sd,log=T)
              mh_factor <- lq_top - lq_bottom
              propose_jpd <- log_posterior(propose_pars,data)
              d <- propose_jpd - current_jpd + sum(mh_factor)
              if (d>=0) {
                  current_pars <- propose_pars
                  current_jpd <- propose_jpd
              } else {
                  if (d>log(runif(1))) {
                      current_pars <- propose_pars
                      current_jpd <- propose_jpd
                  }
              }
              sims.list[[pp]] <- rbind(sims.list[[pp]],c(current_pars[[pp]]))
          }
      }
      ##-#########################################################
      ##   update transition specific IID random effects
      ##-#########################################################
      if (update_setting$update_msm_random) {
          pp <- 'w'
          #   update for each transition
          propose_pars <- current_pars
          lq_top <- lq_bottom <- NULL
          for (jk in ats) {
            which_par <- paste0(pp,'_',jk)
            pms <- paste0(1:fitdata$nctys,'_',jk)
            upt <- gmrf_sampling(which_par,current_pars,data)
            ustar <- upt$x
            propose_pars[[pp]][pms] <- ustar[pms]
            ##  q(ustar|u0)
            tmp <- logden_jump_msm_random(ustar,upt[c('mean','V')])
            lq_bottom <- c(lq_bottom,tmp)
            ##  for q(u0|ustar)
            upt_propose <- gmrf_sampling(which_par,propose_pars,data)
            tmp <- logden_jump_msm_random(current_pars[[pp]][pms],upt_propose[c('mean','V')])
            lq_top <- c(lq_top,tmp)
          }
          data$model_spec$byperson <- F
          propose_jpd <- log_posterior(propose_pars,data)
          d <- propose_jpd - current_jpd + sum(lq_top - lq_bottom)
          if (d>=0) {
              current_pars <- propose_pars
              current_jpd <- propose_jpd
          } else {
              if (d>log(runif(1))) {
                  current_pars <- propose_pars
                  current_jpd <- propose_jpd
              }
          }
          ##   store iteration
          sims.list[[pp]] <- rbind(sims.list[[pp]],c(current_pars[[pp]]))
      }
      ##-#####################################################################
      ##   MH update variance of the transition specific IID random effects
      ##     log(sigma) ~ N(current, jump_sd)
      ##-#####################################################################
      if (update_setting$update_msm_random_SD) {
        jump_sd <- list()
        jump_sd[['sd_w']] <- 0.1
        pp <- 'sd_w'
        propose_pars <- current_pars
        propose_pars[[pp]] <- rlnorm(1,log(current_pars[[pp]]),0.1)
        propose_jpd <- log_posterior(propose_pars,data)
        ##  for q(u_current|u_proposed)
        lq_top <- dlnorm(current_pars[[pp]],propose_pars[[pp]],jump_sd[[pp]],log=T)
        ##  for q(u_proposed|u_current)
        lq_bottom <- dlnorm(propose_pars[[pp]],current_pars[[pp]],jump_sd[[pp]],log=T)
        d <- propose_jpd - current_jpd + lq_top - lq_bottom
        if (d>=0) {
          current_pars <- propose_pars
          current_jpd <- propose_jpd
        } else {
          if (d>log(runif(1))) {
            current_pars <- propose_pars
            current_jpd <- propose_jpd
          }
        }
        sims.list[[pp]] <- c(sims.list[[pp]],current_pars[[pp]])
      }
      
      ##-##################################################
      ##   Gibbs for error SD in LMM
      ##-##################################################
      if (update_setting$update_lmm_data_sd) {
          pp <- 'sd_y'
          ss <- update_data_sd(current_pars,data)
          current_pars[[pp]] <- ss
          sims.list[[pp]] <- c(sims.list[[pp]],ss)
          current_jpd <- log_posterior(current_pars,data)
      }

      ##-##################################################
      ##   Gibbs for hyperparameters in LMM
      ##-##################################################
      if (update_setting$update_lmm_random_effects) {
          ##   update mean vector
          pp <- 'B'
          ss <- sample_lmm_mean_vector(current_pars,data)$x
          current_pars[[pp]] <- matrix(ss,ncol=1)
          sims.list[[pp]] <- rbind(sims.list[[pp]],ss)
          ##   update precision matrix
          pp <- 'Q_D'
          ss <- sample_lmm_precision_matrix(current_pars,data)$x
          current_pars[[pp]] <- ss
          tmp <- array(0,c(1,2,2))
          tmp[1,,] <- ss
          sims.list[[pp]] <- abind::abind(sims.list[[pp]],tmp,along=1)
          ##   calculate the joint posterior at the updated values
          current_jpd <- log_posterior(current_pars,fitdata)
      }
      ##   update progress bar
      if (iter%%(niters/pb_freq)==0) setTxtProgressBar(pb, iter)
  }  ##   next iteration
  ##   close progress bar
  close(pb)
  ##   on screen update
  t2 <- Sys.time()
  ddt <- t2-t1
  msg1 <- paste0(' ** fitting ended at ',t2)
  msg2 <- paste0(' == Fitting time: ',round(ddt,digits=2),' ',attr(ddt,'units'))
  cat(msg1)
  cat(msg2)
  cat(' ')
  return(sims.list)
}

########################################
####  calculate the log joint posterior
########################################
#' @export
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
  prs <- rep(0,20)
  #  age-specific baselines
  prs[1] <- sum(dgamma(params[['l']],shape=0.01,rate=0.01,log=T))
  if (model_spec$weibull_baseline) prs[1] <- prs[1] + sum(dunif(exp(params[['logdel']]),0.00001,5,log=T))  # added on 15Oct2025
  if (model_spec$include_fixed_effects) prs[2] <- sum(dnorm(params[['beta']],mean=0,sd=100,log=T))
  if (model_spec$joint) {
      if (model_spec$lmm_correlated_random_effects) {
          ##   hyperparameters: mean vector
          prs[5] <- emdbook::dmvnorm(c(params[['B']]),rep(0,2),10^2*diag(2),log=T)
          ##   hyperparameters: precision matrix
          rho <- 0
          s1 <- s2 <- 10
          V <- matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2),nrow=2)
          prs[6] <- LaplacesDemon::dwishart(params[['Q_D']],nu=2,S=V,log=T)
          ##   error SD
          prs[7] <- dgamma(1/params[['sd_y']]^2,shape=0.01,rate=0.01,log=T)
          ##  association parameters
          prs[8] <- sum(dnorm(unlist(params[grep('a_',names(params))]),0,100,log=T))
          D <- round(MASS::ginv(params[['Q_D']]),digits=10)
          B <- c(params[['B']])
          b <- cbind(params[['u']],params[['d']])
          ##   correlated intercept-slope
          ##     changed on 08Jul2025 using fastmvn
          ##prs1 <- apply(b,1,function(bi) {dmvnorm(c(bi),B,D,log=T)})  ## before 08Jul2025
          prs1 <- mvnfast::dmvn(b,B,D,log=T)  ## much faster version
          prs2 <- 0
          if (!byperson) prs[9] <- sum(prs1 + prs2)
      } else {
          ##   means of varying intercepts and slopes
          pp <- 'm_u'
          prs[5] <- dnorm(params[[pp]],0,1000,log=T)
          pp <- 'm_d'
          prs[6] <- dnorm(params[[pp]],0,1000,log=T)
          ##   SDs of varying intercepts and slopes
          pp <- 'sd_u'
          prs[7] <- dgamma(1/params[[pp]]^2,shape=0.01,rate=0.01,log=T)
          pp <- 'sd_d'
          prs[8] <- dgamma(1/params[[pp]]^2,shape=0.01,rate=0.01,log=T)
          pp <- 'sd_y'
          prs[9] <- dgamma(1/params[[pp]]^2,shape=0.01,rate=0.01,log=T)
          ##  association parameters
          prs[10] <- sum(dnorm(unlist(params[grep('a_',names(params))]),0,100,log=T))
          if (byperson) {
              ##   random intercepts
              prs1 <- dnorm(params[['u']],mean=params[['m_u']],sd=params[['sd_u']],log=T)
              ##   random slopes
              prs2 <- dnorm(params[['d']],mean=params[['m_d']],sd=params[['sd_d']],log=T)
          } else {
              ##   random intercepts
              prs[11] <- sum(dnorm(params[['u']],mean=params[['m_u']],sd=params[['sd_u']],log=T))
              ##   random slopes
              prs[12] <- sum(dnorm(params[['d']],mean=params[['m_d']],sd=params[['sd_d']],log=T))
          }
      }
  }
  if (model_spec$include_msm_random) {
    #  random effects for each MSM transition
    tmp_prs <- 0
    for (jk in ats) {
      which_pars <- grep(paste0('_',jk),names(params[['w']]))
      pm <- params[['w']][which_pars]
      precision <- 1/params[['sd_w']]^2
      tmp_prs <- tmp_prs + logden_exchangeable_s2z(pm,precision)
    }
    prs[20] <- tmp_prs
    prs[19] <- dunif(params[['sd_w']],0.00001,10,log=T)
  }
  out <- ll + sum(prs)
  if (byperson) out <- out + prs1 + prs2
  return(out)
}
if (F) {
  log_posterior(params,dat)
}


#######################################################
###  function to calculate the fixed effects
###  on the log scale
#######################################################
#' @export
compute_fixed_effects <- function(params, dat, at_quadrature) {
  X <- dat$X
  xvs <- colnames(X)
  ats <- dat$transitions_to_analyse
  trans_indicators <- dat$trans_indicators
  bX <- NULL
  for (k in ats) {
    ##   what betas are involved in this transition
    pbs <- paste0(xvs,'_',k)
    spbs <- intersect(pbs,names(params$beta))  ##  the betas involved in this transition
    xx <- sapply(spbs,function(xxx) {strsplit(xxx,'_')[[1]][1]})
    #  extract beta
    beta <- matrix(params[['beta']][spbs],ncol=1)
    this_X <- X[,xx]
    #  calculate X * beta for this transition
    bX <- cbind(bX,this_X%*%beta*trans_indicators[,k])
##    beta <- matrix(params[['beta']][paste0(dat$XVars,'_',k)],ncol=1)
##    #  calculate X * beta for that transition
##    bX <- cbind(bX,X%*%beta*trans_indicators[,k])
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
#' @export
compute_longitudinal_effects <- function(params, dat, at_quadrature, include_association_par=T) {
    ##   get the relevant data
    dd <- dat
    if (at_quadrature) dd <- dat$Q
    ##pid <- dd$pid
    mergeid <- dd$mergeid
    jk_index <- dd$jk_index
    t_ms_t0 <- dd$t_ms_t0

    ##   extract and format parameters
    ##   a1 and a2 as in a1*CV and a2*CV^2 (CV=current value)
    a1 <- a2 <- rep(0,dat$NTS)
    names(a1) <- names(a2) <- dat$transitions_to_analyse
    a1[names(params[['a_1']])] <- params[['a_1']]
    a2[names(params[['a_2']])] <- params[['a_2']]
    A <- cbind(a1[jk_index],
               a2[jk_index])
    ##   compute (standardised) longitudinal values
    u <- params[['u']][mergeid]
    d <- params[['d']][mergeid]
    sv <- u + d*t_ms_t0
    sv2 <- sv*sv
    ##   CV , CV^2 (CV=current value)
    out <- cbind(sv,sv2)
    if (include_association_par) out <- rowSums(A*out)
    return(out)
}

#######################################################
###  function to calculate the random effects on 
###  MSM transitions on the log scale
#######################################################
#' @export
compute_msm_random_effects <- function(params, dat, at_quadrature) {
    ##   get the relevant data
    dd <- dat
    if (at_quadrature) dd <- dat$Q
    jkc_index <- dd$jkc_index
    out <- params[['w']][jkc_index]
    return(out)
}

#######################################################
###  function to calculate the hazard ratios
#######################################################
#' @export
compute_RR <- function(params, dat, at_quadrature) {
  model_spec <- dat$model_spec
  dd <- dat
  if (at_quadrature) dd <- dat$Q
  bX <- aL <- eta <- w <- 0
  #  fixed effects
  if (model_spec$include_fixed_effects) bX <- compute_fixed_effects(params, dat, at_quadrature)
  #  longitudinal
  if (model_spec$joint) aL <- compute_longitudinal_effects(params, dat, at_quadrature, include_association_par=T)
  #  random effects on MSM transition
  if (model_spec$include_msm_random) w <- compute_msm_random_effects(params, dat, at_quadrature)
  #  putting all together
  logRR <- bX + aL + eta + w
  #  return hazard ratios
  out <- exp(logRR)
  return(out)
}

#######################################################
###  function to calculate the baseline hazards
#######################################################
#' @export
compute_baseline <- function(params, dat, at_quadrature) {
  #  age- and transition-specific baseline hazards
  dd <- dat
  if (at_quadrature) dd <- dat$Q
  h0 <- params[['l']][dd[['jkg_index']]]
  out <- h0
  return(out)
}

#######################################################
###  function to calculate the time-varying part in
###  the Weibull
###               i.e. del*t^(del-1)
###    added on 15Oct2025
#######################################################
#' @export
compute_baseline_weibull_timevarying <- function(params, dat, at_quadrature) {
    #  age- and transition-specific baseline hazards
    dd <- dat
    if (at_quadrature) dd <- dat$Q
    jkg_index <- dd$jkg_index
    t_ms_t0 <- dd$t_ms_t0
    
    del <- exp(params[['logdel']][jkg_index])
    out <- del*t_ms_t0^(del-1)
    return(out)
}

#######################################################
###  function to calculate the transition hazards
###   added h0_tv on 15Oct2025
#######################################################
#' @export
compute_haz <- function(params, dat, at_quadrature, include_quadweights=T, without='') {
  h0 <- compute_baseline(params, dat, at_quadrature)
  h0_tv <- compute_baseline_weibull_timevarying(params, dat, at_quadrature)
  RR <- compute_RR(params, dat, at_quadrature)
  out <- h0*RR*h0_tv
  if (without=='h0') out <- RR*h0_tv
  if (at_quadrature & include_quadweights) out <- out * dat$Q$w  #  multiplying quadrature weights
  return(out)
}

#######################################################
###  function to calculate the multistate likelihood
#######################################################
#' @export
ll_msm <- function(params, dat, byperson=FALSE) {
  status <- dat[['status']]
  ll1 <- compute_haz(params,dat,at_quadrature = F)^status
  ll2 <- compute_haz(params,dat,at_quadrature = T)  # this includes quadrature weights
  if (byperson) {
    tll1 <- tapply(log(ll1),dat$mergeid,sum)
##    ids <- sapply(dat$subject,function(x){which(names(tll1)==x)})
##    tll1 <- tll1[ids]
    tll1 <- tll1[dat$subject]
    tll2 <- tapply(ll2,dat$Q$mergeid,sum)
##    ids <- sapply(dat$subject,function(x){which(names(tll2)==x)})
##    tll2 <- tll2[ids]
    tll2 <- tll2[dat$subject]
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
#' @export
ll_lmm <- function(params, dat, byperson=F) {
  dd <- dat$long
  y <- dd$y
  tt <- dd$t_ms_t0
##  pid <- dd$pid
  mergeid <- dd$mergeid
  u <- params[['u']][mergeid] # intercepts
  d <- params[['d']][mergeid] # slopes
  m <- u + d*tt
  sy <- params[['sd_y']]
  ll <- dnorm(y,mean=m,sd=sy,log=T)
  if (byperson) {
    out <- tapply(ll,mergeid,sum)
    ##ids <- sapply(dat$subject,function(x){which(names(out)==x)})
    ##out <- out[ids]
    out <- out[dat$subject]
  } else {
    out <- sum(ll)
  }
  return(out)
}

#######################################################
###  function to update age-specific baselines
#######################################################
#' @export
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
###  function to update transition specific random effects
#######################################################
#' @export
update_msm_random_effects <- function(params, dat) {
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
###  identify the function to approximate the density
#######################################################
#' @export
identify_approx_function <- function(which_par) {
    out <- NULL
    if (length(grep('a_',which_par))>0) out <- f_association_approximated
    if (which_par=='beta') out <- f_beta_approximated
    if (which_par=='intercept' | which_par=='slope') out <- f_lmm_random_effects_approximated
    if (which_par=='eta') out <- f_eta_approximated
    if (length(grep('w_',which_par))>0) out <- f_w_approximated
    ##   added 12Jun2025
    if (which_par=='lmm_corr_random_effects') out <- f_lmm_corr_random_effects_approximated
    if (is.null(out)) stop('no approximation function has been selected')
    return(out)
}
    
#######################################################
###  function to find the mode iteratively
#######################################################
#' @export
auto_mode_finder <- function(which_par, start_x0, params, dat) {
    ##
    ##   updated on 31Jan2025: fixed the problem that v2_31Jan2025 returns
    ##                         a value found at the penultimate step, not the
    ##                         final optimal one
    ##   updated on 05Jun2025: using Newton-Raphson to find the maximum of
    ##                         each density (see p. 169 Rue and Held 2005)
    ##                         there is no need to compute the density
    ##
    gmrf_settings <- dat$gmrf_settings
    xr <- gmrf_settings$search_range
    sp <- gmrf_settings$search_points
    xv <- seq(xr[1],xr[2],length.out=sp)
    
    ##   identify the function to approximate the density
    f <- identify_approx_function(which_par)

    ##   initialise counter
      iter <- 0
      ##   get starting values
      x0 <- NULL
      if (which_par!='lmm_corr_random_effects') {
          if (length(grep('w_',which_par))>0) {
            x0 <- start_x0
          } else {
            nm <- names(start_x0)
            x0 <- start_x0
          }
      } else {
          sb <- dat$subject
          x0 <- cbind(params[['u']][sb],params[['d']][sb])
          rownames(x0) <- sb
      }
      current_x0 <- x0
      continue <- TRUE
      store_x0 <- NULL
      ##   carry out Newton-Raphson
      while (continue) {
        ##   get the new positions
          ##--- before 05Jun2025
        ## d <- f(xv,x0,params,dat,which_par)$d
        ## new_x0 <- apply(d,2,function(xx){xv[which.max(xx)]})
        ##   modified on 05Jun2025
        ff <- f(xv,current_x0,params,dat,which_par)
        new_x0 <- ff$eta_next
        if (which_par!='lmm_corr_random_effects') {
            if (length(grep('w_',which_par))>0) {
                m <- ff$m
                V <- ff$V
            } else {
                names(new_x0) <- nm
                dev1 <- ff$dev1
                dev2 <- ff$dev2
            }
        } else {
            m <- ff$m
            V <- ff$V
        }
        ##   calculate absolute relative error (in 100%)
        err <- c(abs((current_x0-new_x0)/new_x0)*100)
        ##   update counter
        iter <- iter + 1
        ##   continue or stop
        if (all(err[which(new_x0!=0)]<=gmrf_settings$eps)) {
            ##   all reached maximum
            continue <- FALSE
        } else {
            ##   some have not reached maximum
            if (iter==gmrf_settings$max_iters) {
                ##   max. iteration reached
                continue <- FALSE
                stop('max iteration reached but some still not at maximum')
            }
        }
        ##   continue or not, update the positions
        current_x0 <- new_x0
        store_x0 <- rbind(store_x0,new_x0)
      }
    if (which_par!='lmm_corr_random_effects') {
      if (length(grep('w_',which_par))>0) {
          out <- list(x=current_x0,iter=iter,abs_err100=err,m=m,V=V)
      } else {
          out <- list(x=current_x0,iter=iter,abs_err100=err,dev1=dev1,dev2=dev2)
      }
    } else {
        out <- list(x=current_x0,iter=iter,abs_err100=err,m=m,V=V)
    }
    return(out)
}

#' @export
sample_at_opt_mode_univariate <- function(opt) {
  mx <- opt$x
  dev1 <- opt$dev1
  dev2 <- opt$dev2
  out <- approx_coefficients(NULL,mx,dev1,dev2,calc_den=F)
  out$x <- c(sapply(1:length(mx),function(xx) {rnorm(1,out$mean[xx],out$sd[xx])}))
  return(out)
}

#' @export
sample_at_opt_mode_multivariate <- function(opt) {
  mx <- opt$x
  V <- opt$V
  s <- c(MASS::mvrnorm(1,mx,V))
  out <- list(mean=mx,V=V,x=s)
  return(out)
}

#' @export
check_positive_definiteness <- function(V) {
  if (length(dim(V))==2) {
    tmp <- array(0,c(1,dim(V)))
    tmp[1,,] <- V
    V <- tmp
  }
  c1 <- which(signif(V[,1,2],digits=10)!=signif(V[,2,1],digits=10))
  dtm <- V[,1,1]*V[,2,2] - V[,1,2]*V[,2,1]
  c2 <- which(dtm<=0)
  c3 <- which(V[,1,1]<=0)
  c4 <- which(V[,2,2]<=0)
  msg <- NULL
  out <- list()
  ie <- 0
  if (length(c1)>0) {
    ie <- ie + 1
    msg <- c(msg,'some V are not symmetric')
    out[[ie]] <- c1
  }
  if (length(c2)>0) {
    ie <- ie + 1
    msg <- c(msg,'some V have non-positive determinant')
    out[[ie]] <- c2
  }
  if (length(c3)>0) {
    ie <- ie + 1
    msg <- c(msg,'some V have non-positive variance for first component')
    out[[ie]] <- c3
  }
  if (length(c4)>0) {
    ie <- ie + 1
    msg <- c(msg,'some V have non-positive variance for second component')
    out[[ie]] <- c4
  }
  out <- list(msg=msg,out=out)
  return(out)
}

#' @export
sample_at_opt_mode_multivariate_faster <- function(opt) {
  ##   https://www.probabilitycourse.com/chapter5/5_3_2_bivariate_normal_dist.php
  mx <- opt$x
  V <- opt$V
  m_x <- opt$x[,1]  ##   the means
  m_y <- opt$x[,2]  ##   the means
  sd_x <- sqrt(V[,1,1])
  sd_y <- sqrt(V[,2,2])
  rhos <- V[,1,2]/(sd_x*sd_y)  ##   corr. coeffs
  n <- length(m_x)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  sx <- cbind(sd_x*z1 + m_x,
              sd_y*(rhos*z1+sqrt(1-rhos^2)*z2) + m_y)
  V <- opt$V
    out <- list(mean=mx,V=V,x=sx)
  return(out)
}


#######################################################
###  function to sample parameter values based on GMRF
#######################################################
#' @export
gmrf_sampling <- function(which_par, params, dat) {
    if (which_par=='intercept') pp <- 'u'
    if (which_par=='slope') pp <- 'd'
    if (length(grep('a_',which_par))>0) pp <- which_par
    if (which_par=='beta') pp <- which_par
    is_cty_random <- FALSE
    if (length(grep('w_',which_par))>0) is_cty_random <- TRUE
    if (which_par=='lmm_corr_random_effects') {
        pp <- which_par
        start_x0 <- NULL
    } else {
        if (is_cty_random) {
            #  random effects on an MSM transition
            pp <- 'w'
            which_jk <- sub(pp,'',which_par)
            start_x0 <- params[[pp]][grep(which_jk,names(params[['w']]))]
        } else {
            start_x0 <- params[[pp]]
        }
    }

    ##   find the modes via Newton-Raphson
    opt_x <- auto_mode_finder(which_par, start_x0, params, dat)
    ##   sample based on the modes
    if (which_par!='lmm_corr_random_effects') {
      if (!is_cty_random) {
        xs <- sample_at_opt_mode_univariate(opt_x)
        names(xs$x) <- names(start_x0)
        out <- list(mean=xs$mean,sd=xs$sd,x=xs$x)
      } else {
        v_check <- check_positive_definiteness(opt_x$V)  ##  check positive definiteness
        if (!is.null(v_check$msg)) report_problems(v_check, dat)
        xs <- sample_at_opt_mode_multivariate(opt_x)
        out <- list(mean=xs$mean,V=xs$V,x=xs$x)
      }
    } else {
      v_check <- check_positive_definiteness(opt_x$V)  ##  check positive definiteness
      if (!is.null(v_check$msg)) report_problems(v_check, dat)
      xs <- sample_at_opt_mode_multivariate_faster(opt_x)
      out <- list(mean=xs$mean,V=xs$V,x=xs$x)
    }
    return(out)
}

#' @export
report_problems <- function(vc, dat) {
  for (i in 1:length(vc$msg)) {
    print(vc$msg[i])
    print(vc$out[[i]])
  }
  pbs <- unique(unlist(lapply(vc$out,names)))
  for (ip in pbs) {
    print(paste0('--- ',ip,' ---'))
    ##   msm data
    print(' MSM data')
    print(as.data.frame(dat[c('mergeid','from','to','start','stop','status','age_entry')])[which(dat$mergeid==ip),])
    print(' long data')
    print(as.data.frame(dat$long[c('mergeid','y','age','t_ms_t0')])[which(dat$long$mergeid==ip),])
  }
  stop(' problem with some of the covariance matrices')
}

#######################################################
###  function to calculate the coefficients in the
###  approximation aa+bb*x-1/2cc*x^2
###  with aa, bb and cc calculated at x0
###
###  aa <- f0 - x0*dev1 + x0^2/2*dev2
#######################################################
#' @export
approx_coefficients <- function(x,x0,dev1,dev2,calc_den=F) {
    ##   modified on 05Jun2025: added calc_den to suppress the calculation
    ##                          of the density at a sequence of values
    ##                          that's not needed because the maximum is found
    ##                          via Newton-Raphson
    bb <- dev1 - x0*dev2
    cc <- -dev2
    mm <- bb/cc  # mean
    prec <- cc   # precision
    ################################################
    ##    the whole calc_den part is redundant!!
    ################################################
    d <- 0
    if (calc_den) {
        ##   approximated density at x using Taylor expansion at x0
        d <- matrix(0,nrow=length(x),ncol=length(x0))
        ##for (ix in 1:length(x0)) d[,ix] <- f0[ix] + dev1[ix]*(x-x0[ix]) + dev2[ix]/2*(x-x0[ix])^2
        for (ix in 1:length(x0)) d[,ix] <- dev1[ix]*(x-x0[ix]) + dev2[ix]/2*(x-x0[ix])^2
    }
    ##   putting everything together for output
    out <- list(b=bb,c=cc,d=0)
    out$mean <- mm
    out$sd <- sqrt(1/prec)
    return(out)
}

#' @export
update_random_effect_mean <- function(params,dat,which_lmm_part) {
  #
  #  prior is in the form of N(prior[1]=0,prior[2]=precision)
  #
  #   which_random_effect <- 'u'  #  intercept
  #   which_random_effect <- 'd'  #  slope
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

#' @export
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

#' @export
update_data_sd <- function(params,dat) {
  prior <- c(0.01,0.01)
  dl <- dat$long
  y <- dl$y      #  longitudinal observations
  tt <- dl$t_ms_t0
  nts <- dl$nts  #  number of observations per person
##  pid <- dl$pid  #  which y belongs to which person
  u <- params[['u']][dl$mergeid] # intercept
  d <- params[['d']][dl$mergeid] # slope
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


#' @export
sample_lmm_mean_vector <- function(params, dat) {
    p <- 2  #  number of correlated random effects per person
    npds <- dat$npds
    u <- params[['u']]
    d <- params[['d']]
    ##   mean vector of the current random effects
    b <- matrix(c(mean(u),mean(d)),nrow=p)
    Q_D <- params[['Q_D']]
    ##   prior specification
    B_0 <- matrix(rep(0,p),nrow=p)
    Q_D0 <- 1/100^2*diag(p)
    ##   mean and variance of the full conditional
    P <- npds*Q_D + Q_D0
    V <- MASS::ginv(P)
    M <- t(npds*Q_D%*%b + Q_D0%*%B_0)%*%V
    x <- MASS::mvrnorm(1,mu=M,Sigma=V)
    out <- list(M=M,V=V,x=x)
    return(out)
}

#' @export
sample_lmm_precision_matrix <- function(params, dat) {
    ##   Wishart prior
    rho <- 0
    s1 <- s2 <- 10
    V <- matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2),nrow=2)
    inv_V <- MASS::ginv(V)
    nu0 <- 2

    npds <- dat$npds
    u <- params[['u']]
    d <- params[['d']]
    b <- cbind(u,d)
    B <- params[['B']]
    db <- b
    for (i in 1:ncol(b)) db[,i] <- b[,i] - B[i]

    S <- MASS::ginv(t(db)%*%db + inv_V)
    nu <- npds + nu0
    x <- LaplacesDemon::rwishart(nu, round(S,digits=11))
    out <- list(roundS=round(S,digits=10),S=S,nu=nu,x=x)
    return(out)
}
