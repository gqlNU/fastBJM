rm(list=ls())

library(statmod)
workdir <- '/Volumes/WorkSpace/SHARE/data/health/'
setwd(workdir)
source('../../code/health/functions.R')
source('../../code/health/functions4fits.R')

source('/Volumes/WorkSpace/OnGitHub/fastBJM/functions.r')
source('/Volumes/WorkSpace/OnGitHub/fastBJM/functions4fits.r')
source('/Volumes/WorkSpace/OnGitHub/fastBJM/functions4approximation.r')
source('/Volumes/WorkSpace/OnGitHub/fastBJM/functions4IGMRF.r')
############################
##   user input
############################
niters <- 100
nrefresh <- 1
outcome <- 'hypertension_stroke'
data_type <- 'joint'
icty <- c(1:8,13)
icty <- 1:3
gd <- 1                 #  gender 1=male 2=female
X_spec <- matrix(c('marital','cat',
                   'education','cnt'),byrow=T,ncol=2)
X_spec <- matrix(c('education','cnt'),byrow=T,ncol=2)
nQs <- 1  #  number of quadrature points
age_range <- c(50,150)  #  lower and upper bounds in defining the baseline age intervals
age_cuts <- c(60,70,80)
allowable_transitions <- ats <- c('12','13','15','24','25','34','35','45')
transitions_to_analyse <- c('12')
which_longvar <- 'bmi'

model_spec <- list()
model_spec$msm_only <- F        #  only modelling the multistate outcome if T; else joint modelling
model_spec$msm_by_country <- ''
model_spec$msm_country_age <- FALSE
##   end user input
############################
model_spec <- gather_model_spec(model_spec, X_spec)
update_setting <- gather_update_setting(model_spec)

fitdata <- prepare_fitdata(icty, gd, data_type, 
                           allowable_transitions, transitions_to_analyse,
                           X_spec,
                           age_range, age_cuts,
                           nQs,
                           which_longvar)
fitdata$model_spec <- model_spec

params <- get_parameters(fitdata, model_spec)

params[['l']] <- update_baseline(params, fitdata, c(0.01,0.01))
if (update_setting$update_fixed_effects) params[['beta']] <- auto_mode_finder('beta', params[['beta']], params, fitdata)$x0
if (update_setting$update_lmm_random_effects) {
  params[['u']] <- auto_mode_finder('intercept', params[['u']], params, fitdata)$x0
  params[['d']] <- auto_mode_finder('slope', params[['d']], params, fitdata)$x0
}
if (update_setting$update_association) {
  for (pp in c('a_1','a_2')) params[[pp]] <- auto_mode_finder(pp, params[[pp]], params, fitdata)$x0
}
if (update_setting$update_cty_random_effects) params[['eta']] <- sample_country_random_effects(params, fitdata, eval=NULL)$Ax


inits1 <- params
t1 <- tstart <- Sys.time()
store <- NULL
npars <- length(inits1)
current_pars <- inits1
current_jpd <- log_posterior(current_pars, fitdata)

sims.list <- as.list(1:npars)
names(sims.list) <- names(inits1)
for (i in 1:npars) {
  sims.list[[i]] <- current_pars[[i]]
}

print(fitdata$ncases_by_jkg)

model_spec$byperson <- F
ats <- fitdata$transitions_to_analyse




for (iter in 1:niters) {
  ####################################################
  #  age-specific baselines
  ####################################################
  if (update_setting$update_h0) {
    pp <- 'l'
    x <- update_baseline(current_pars, fitdata, prior=c(0.01,0.01))
    nms <- names(x)
    current_pars[[pp]][nms] <- x[nms]
    sims.list[[pp]] <- rbind(sims.list[[pp]],current_pars[[pp]])
    current_jpd <- log_posterior(current_pars, fitdata)
  }
  
  ###########################################################
  #   update intercepts and slopes (p. 169 in Rue and Held)
  ###########################################################
  if (update_setting$update_lmm_random_effects) {
    for (which_lmm_part in c('intercept','slope')) {
      model_spec$byperson <- T
      propose_pars <- current_pars
      current_jpd_byperson <- log_posterior(current_pars, fitdata)
      if (which_lmm_part=='intercept') pp <- 'u'
      if (which_lmm_part=='slope') pp <- 'd'
      u0 <- current_pars[[pp]]				# current u
      upt <- gmrf_sampling(which_lmm_part,current_pars,fitdata)
      ustar <- upt$x
      propose_pars[[pp]] <- ustar  	# proposed u
      propose_jpd_byperson <- log_posterior(propose_pars,fitdata)
      #  q(ustar|u0)
      ind_lq_bottom <- sapply(1:fitdata$npds,function(xx){dnorm(ustar[xx],mean=upt$mean[xx],sd=upt$sd[xx],log=T)})
      lq_bottom <- sum(ind_lq_bottom)
      #  for q(u0|ustar)
      upt_propose <- gmrf_sampling(which_lmm_part,propose_pars,fitdata)      
      ind_lq_top <- sapply(1:fitdata$npds,function(xx){dnorm(u0[xx],mean=upt_propose$mean[xx],sd=upt_propose$sd[xx],log=T)})
      lq_top <- sum(ind_lq_top)
      d <- propose_jpd_byperson - current_jpd_byperson + ind_lq_top - ind_lq_bottom
      lr <- log(runif(length(d)))
      icc <- which(d>=lr)
      current_pars[[pp]][icc] <- propose_pars[[pp]][icc]
      model_spec$byperson <- F
      current_jpd <- log_posterior(current_pars,fitdata)
      sims.list[[pp]] <- rbind(sims.list[[pp]],current_pars[[pp]])
    }
  }
  model_spec$byperson <- F

  ############################
  #  update fixed effects
  ############################
  if (update_setting$update_fixed_effects) {
    pp <- 'beta'
    propose_pars <- current_pars
    bstar <- gmrf_sampling(pp, current_pars, fitdata)
    propose_pars[[pp]] <- bstar$x
    #  q(ustar|u0)
    lq_bottom <- sum(dnorm(bstar$x,mean=bstar$mean,sd=bstar$sd,log=T))
    bc <- gmrf_sampling(pp, propose_pars, fitdata)
    #  q(u0|ustar)
    lq_top <- sum(dnorm(current_pars[[pp]],mean=bc$mean,sd=bc$sd,log=T))
    mh_factor <- lq_top - lq_bottom
    propose_jpd <- log_posterior(propose_pars, fitdata)
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
  
  ##################################
  #  update association parameters
  ##################################
  if (update_setting$update_association) {
    for (pp in c('a_1','a_2')) {
      propose_pars <- current_pars
      bstar <- gmrf_sampling(pp, current_pars, fitdata)
      propose_pars[[pp]] <- bstar$x
      #  q(ustar|u0)
      lq_bottom <- dnorm(bstar$x,mean=bstar$mean,sd=bstar$sd,log=T)
      bc <- gmrf_sampling(pp, propose_pars, fitdata)
      #  q(u0|ustar)
      lq_top <- dnorm(current_pars[[pp]],mean=bc$mean,sd=bc$sd,log=T)
      mh_factor <- lq_top - lq_bottom
      propose_jpd <- log_posterior(propose_pars,fitdata)
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

  #########################################
  #  update country-level random effects
  #########################################
  if (update_setting$update_cty_random_effects) {
    propose_pars <- current_pars
    pp <- 'eta'
    v0 <- current_pars[[pp]]				# current eta
    vpt <- sample_country_random_effects(current_pars, fitdata,eval=NULL)
    vstar <- vpt$Ax
    uncons <- vpt$x
    propose_pars[[pp]] <- vstar  	# proposed eta
    propose_jpd <- log_posterior(propose_pars,fitdata)
    #  q(vstar|v0)
    lq_bottom <- vpt$log_den
    #  for q(u0|ustar)
    vct <- sample_country_random_effects(propose_pars,fitdata,eval=list(x=uncons,Ax=v0))
    lq_top <- vct$log_den
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
    sims.list[[pp]] <- rbind(sims.list[[pp]],c(current_pars[[pp]]))
  }
  ####################################################
  #   MH for country random effect SDs
  #    option 1: Eq. 4.13 on p. 146 in Rue and Held
  #    option 2: looking into the full conditionals?
  ####################################################
  if (update_setting$update_cty_random_effects_SD) {
    pp <- 'sd_eta'
    for (k in ats) {
      propose_pars <- current_pars
      propose_pars[[pp]][k] <- runif(1,0.01,2)  #  the lower bound cannot be 0 to guard against the problem of 
      #  the values in the first eigenvector being unequal
      #  to see that m <- array(1e-10,c(20,20));diag(m) <- 1/0.001^2;eigen(m)$vectors[,1]
      #  then IGMRF_Q will return stop('something strange in the eigenvectors')
      propose_jpd <- log_posterior(propose_pars,fitdata)
      d <- propose_jpd - current_jpd
      if (d>=0) {
        current_pars <- propose_pars
        current_jpd <- propose_jpd
      } else {
        if (d>log(runif(1))) {
          current_pars <- propose_pars
          current_jpd <- propose_jpd
        }
      }
    }
    sims.list[[pp]] <- rbind(sims.list[[pp]],current_pars[[pp]])
  }
  
  #########################################
  #  update age-country interactions
  #########################################
  if (update_setting$update_age_cty_interactions) {
    propose_pars <- current_pars
    pp <- 'kappa'
    v0 <- current_pars[[pp]]				# current kappa
    vpt <- sample_age_country_interactions(current_pars,fitdata,eval=NULL)
    vstar <- vpt$Ax
    uncons <- vpt$x
    propose_pars[[pp]] <- vstar  	# proposed kappa
    propose_jpd <- log_posterior(propose_pars,fitdata)
    #  q(vstar|v0)
    lq_bottom <- vpt$log_den
    #  for q(u0|ustar)
    vct <- sample_age_country_interactions(propose_pars,fitdata,eval=list(x=uncons,Ax=v0))
    lq_top <- vct$log_den
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
    sims.list[[pp]] <- rbind(sims.list[[pp]],c(current_pars[[pp]]))
  }
  ####################################################
  #   MH for SD of age-country interactions
  #    option 1: Eq. 4.13 on p. 146 in Rue and Held
  #    option 2: looking into the full conditionals?
  ####################################################
  if (update_setting$update_age_cty_interactions_SD) {
    pp <- 'sd_kappa'
    propose_pars <- current_pars
    propose_pars[[pp]] <- runif(1,0.01,2)  #  the lower bound cannot be 0 to guard against the problem of 
    #  the values in the first eigenvector being unequal
    #  to see that m <- array(1e-10,c(20,20));diag(m) <- 1/0.001^2;eigen(m)$vectors[,1]
    #  then IGMRF_Q will return stop('something strange in the eigenvectors')
    propose_jpd <- log_posterior(propose_pars,fitdata)
    d <- propose_jpd - current_jpd
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
  
  ####################################################
  #   Gibbs for mu_u, s and su in LMM
  ####################################################
  if (update_setting$update_lmm_random_effects | update_setting$update_lmm_data_sd) {
    for (pp in c('m_u','sd_u','m_d','sd_d','sd_y')) {
      if (update_setting$update_lmm_random_effects) {
        if (pp=='m_u')  ss <- update_random_effect_mean(current_pars,fitdata,'intercept')
        if (pp=='m_d')  ss <- update_random_effect_mean(current_pars,fitdata,'slope')
        if (pp=='sd_u') ss <- update_random_effect_sd(current_pars,fitdata,'intercept')
        if (pp=='sd_d') ss <- update_random_effect_sd(current_pars,fitdata,'slope')
      }
      if (update_setting$update_lmm_data_sd) {
        if (pp=='sd_y') ss <- update_data_sd(current_pars,fitdata)
      }
      current_pars[[pp]] <- ss
      sims.list[[pp]] <- c(sims.list[[pp]],ss)
    }
    current_jpd <- log_posterior(current_pars,fitdata)
  }
  if (iter%%nrefresh==0) print(paste0(iter,' out of ',niters))
}