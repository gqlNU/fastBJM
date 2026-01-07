#' @export
gather_update_setting <- function(dat, model_spec) {
  out <- list()
  ##   update baseline hazards?
  out$update_h0 <- T
  ##   update fixed effect(s)?
  out$include_fixed_effects <- model_spec$include_fixed_effects
  out$update_fixed_effects <- out$include_fixed_effects
  ##   update association parameter(s)?
  out$update_association <- T
  out$update_lmm_random_effects <- T
  out$update_lmm_data_sd <- T
  if (!model_spec$joint) {
    out$update_association <- F
    out$update_lmm_random_effects <- F
    out$update_lmm_data_sd <- F
  }
  ##   which fixed effects to be updated during MCMC
  bb <- NULL
  msb <- model_spec$beta_by_country
  for (ix in 1:length(model_spec$beta_by_country)) {
    xn <- names(msb)[ix]
    xtype <- dat$X_spec$type[which(dat$X_spec$xn==xn)]
    if (xtype=='cnt') {
        ##   a continuoue-valued covariate
        tmp <- paste0(xn,'_',msb[[ix]])
    } else {
        ##   a categorical covariate
        tmp <- sapply(msb[[ix]],function(xx) {paste0(dat$XVars[grep(xn,dat$XVars)],'_',xx)})
    }
    bb <- c(bb,tmp)
  }
  out$beta_to_update <- bb

  ##   update weibull scale?
  out$update_logdel <- FALSE
  if (model_spec$weibull_baseline) out$update_logdel <- TRUE
  return(out)
}



########################################
####  prepare parameters
########################################
#' @export
get_parameters <- function(dat, model_spec) {
    ats <- dat$transitions_to_analyse
    ps <- c(sapply(ats,function(k){paste0(k,'g',1:dat$nGs)}))
    ##########################
    ##   baseline
    ##########################
    ##   added on 15Oct2025
    params <- list(l=rep(0.05,length(ps)),logdel=rep(0,length(ps)))
    names(params[['l']]) <- names(params[['logdel']]) <- ps
    
    ##########################
    ##   fixed effects
    ##########################
    if (model_spec$include_fixed_effects) {
        bc <- model_spec$beta_by_country
        nms <- NULL
        for (ibc in 1:length(bc)) {
            xn <- names(bc)[ibc]
            xtype <- dat$X_spec$type[which(dat$X_spec$xn==xn)]
            if (xtype=='cnt') {
                ##   a continous-valued covariate
                tmp <- paste0(xn,'_',bc[[ibc]])
            } else {
                ##   a categorical covariate
                tmp <- sapply(bc[[ibc]],function(xx) {paste0(dat$XVars[grep(xn,dat$XVars)],'_',xx)})
            }
            nms <- c(nms,tmp)
        }
        params[['beta']] <- rep(0,length(nms))
        names(params[['beta']]) <- nms
    }
    if (model_spec$joint) {
        ##########################
        ##   association parameters
        ##########################
        params[['a_1']] <- rep(0,length(model_spec$a1p))
        names(params[['a_1']]) <- model_spec$a1p
        params[['a_2']] <- rep(0,length(model_spec$a2p))
        names(params[['a_2']]) <- model_spec$a2p
        ##   error SD
        params[['sd_y']] <- sd(dat$long$y)
        ##   random intercepts and slopes
        params[['u']] <- rep(mean(dat$long$y),dat$npds)
        params[['d']] <- rep(0,dat$npds)
        names(params$u) <- names(params$d) <- dat$subject
        ##   hyperparameters
        if (model_spec$lmm_correlated_random_effects) {
            ##   mean vector and precision matrix
            params[['B']] <- matrix(c(mean(dat$long$y),0),nrow=2)
            params[['Q_D']] <- diag(2)
        } else {
            ##   independent random effects
            params[['sd_u']] <- 1
            params[['sd_d']] <- 1
            params[['m_u']] <- 0
            params[['m_d']] <- 0
        }
    }
    return(params)
}



########################################
####  prepare parameters
########################################
#' @export
initialise_parameters <- function(dat, params, model_spec, update_setting) {
    out <- params
    ##   correlated random effects
    if (update_setting$update_lmm_random_effects) {
        pp <- 'lmm_corr_random_effects'
        tmp <- gmrf_sampling(pp, params, dat)
        out[['u']] <- tmp$x[,1]
        out[['d']] <- tmp$x[,2]
        names(out[['u']]) <- names(out[['d']]) <- dat$subject
        out[['B']] <- sample_lmm_mean_vector(out,dat)$x
        out[['Q_D']] <- sample_lmm_precision_matrix(out,dat)$x
        out[['sd_y']] <- update_data_sd(out,dat)
    }
    out[['l']] <- update_baseline(out, dat, c(0.01,0.01))
    if (update_setting$update_fixed_effects)
        out[['beta']] <- gmrf_sampling('beta', out, dat)$x
    return(out)
}


########################################
####  setup storage for MCMC iterations
########################################
#' @export
setup_storage <- function(params) {
    store <- NULL
    npars <- length(params)
    store <- as.list(1:npars)
    names(store) <- names(params)
    for (i in 1:npars) {
        if (names(store)[i]=='Q_D') {
            tmp <- array(0,c(1,2,2))
            tmp[1,,] <- params[[i]]
            store[[i]] <- tmp
        } else if (names(store)[i]=='B') {
            store[[i]] <- c(params[[i]])
        } else {
            store[[i]] <- params[[i]]
        }
    }
    return(store)
}

