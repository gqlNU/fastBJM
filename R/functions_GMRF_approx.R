########################################
####  fixed effects
########################################
#' @export
extract_jk_x <- function(xn) {
	tmp <- strsplit(xn,'_')[[1]]
	out <- list(jk=tmp[2],xc=tmp[1])
	return(out)
}


#' @export
f_beta_approximated <- function(bx,bx0,params,dat,which_par='') {
  QD <- dat$Q
  pm0 <- params
  nms <- names(bx0)
  pm0[['beta']][nms] <- bx0[nms]
  status <- dat$status
  dev1 <- dev2 <- rep(0,length(bx0))
  names(dev1) <- names(dev2) <- nms

  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  hr <- as.numeric(tapply(hq,dat$Q$row_id,sum))
  ##   go through each beta
  for (xn in nms) {
  	this_x <- extract_jk_x(xn)  ##  extract name and jk
  	jk <- this_x$jk
  	xc <- this_x$xc
  	tv <- dat$trans_indicators[,jk]  ##  which transition it's involved with
  	pp <- xn
    xt <- dat$X[,xc]  ##   the X values
    #  first and second derivatives of the full conditionals calculated at bx0
    dvll <- sum(status*tv*xt) - sum(hr*tv*xt)
    dev1[pp] <- dvll - 1/100^2*bx0[pp]
    dvll <- - sum(hr*tv*xt^2)
    dev2[pp] <- dvll - 1/100^2
  }
  ##   added on 05Jun2025
  ##   Newton-Raphson to find eta such that
  ##   the first derivative is 0 (see p. 169 Rue and Held 2005)
  eta_next <- bx0 - dev1/dev2
  out <- list(eta_next=eta_next,dev1=dev1,dev2=dev2)
  return(out)
}

########################################
####  random effects in longitudinal
########################################
#' @export
f_lmm_random_effects_exact <- function(bx,bx0,params,dat,which_lmm_part) {
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

    ##   from MSM
    ##    hazard at quadrature points
    hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
    ##    longitudinal measurements at event times
    lv <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,1]
    ##    longitudinal measurements at quadrature
    lvq <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,1]

	full_a1 <- full_a2 <- rep(0,dat$NTS)
    names(full_a1) <- names(full_a2) <- dat$transitions_to_analyse
    full_a1[names(pm0[['a_1']])] <- pm0[['a_1']]
    full_a2[names(pm0[['a_2']])] <- pm0[['a_2']]
    
    a1 <- full_a1[dat$jk_index]
    a1q <- full_a1[dat$Q$jk_index]
    a2 <- full_a2[dat$jk_index]
    a2q <- full_a2[dat$Q$jk_index]

  #  f from MSM
  #    the quadrature part
  ##ahq <- (a1q*t_ms_t0q + 2*a2q*t_ms_t0q*lvq) * hq
  ahq <- hq
  ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
  #    the event time part
  ##sat <- status*(a1*t_ms_t0 + 2*a2*t_ms_t0*lv)
  sat <- status*(a1*lv + a2*lv^2)
  sat_p <- tapply(sat,dat$pid,sum)
  #    putting the two together
  dev0_msm <- sat_p - ahq_p

    ##   from longitudinal
    dl <- dat$long
    y <- dl$y
    tau_y <- 1/params[['sd_y']]^2
    tt <- rep(1,length(dl$pid))
    if (pp=='d') tt <- dl$t_ms_t0

    u <- pm0[['u']][dl$pid]
    d <- pm0[['d']][dl$pid]
    mm <- u + d*dl$t_ms_t0
    sd_y <- 1/sqrt(tau_y)
    d0 <- sapply(1:length(y),function(ii){dnorm(y[ii],mm[ii],sd_y,log=T)})
    dev0_long <- as.numeric(tapply(d0,dl$pid,sum))

  ###  from prior
  tau_RE <- 1/pm0[[paste0('sd_',pp)]]^2
  sd_RE <- pm0[[paste0('sd_',pp)]]
  mean_RE <- pm0[[paste0('m_',pp)]]
  RE <- pm0[[pp]]
  dev0_prior <- sum(dnorm(RE,mean_RE,sd_RE,log=T))

  ###  putting both parts together
  f0 <- dev0_msm + dev0_long + dev0_prior
  out <- f0
  return(out)
}



########################################
####  random effects in longitudinal
########################################
#' @export
f_lmm_random_effects_approximated <- function(bx,bx0,params,dat,which_lmm_part) {
	##   ===  updated on 05Jun2025
	##   rectified the mistake of not carrying out the product rule
	##   in deriving the second derivative
	##   (i.e. d(f(d[i])*exp(g(d[i])))/d(d[i]))
	##
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

    ##   from MSM
    ##    hazard at quadrature points
    hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
    ##    longitudinal measurements at event times
    lv <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,1]
    ##    longitudinal measurements at quadrature
    lvq <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,1]

	full_a1 <- full_a2 <- rep(0,dat$NTS)
    names(full_a1) <- names(full_a2) <- dat$transitions_to_analyse
    full_a1[names(pm0[['a_1']])] <- pm0[['a_1']]
    full_a2[names(pm0[['a_2']])] <- pm0[['a_2']]
    
    a1 <- full_a1[dat$jk_index]
    a1q <- full_a1[dat$Q$jk_index]
    a2 <- full_a2[dat$jk_index]
    a2q <- full_a2[dat$Q$jk_index]
    
##    a1 <- pm0[['a_1']][dat$jk_index]
##    a1q <- pm0[['a_1']][dat$Q$jk_index]
##    a2 <- pm0[['a_2']][dat$jk_index]
##    a2q <- pm0[['a_2']][dat$Q$jk_index]

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
  ##  incorrect version:  ahq <- (2*a2q*t_ms_t02q) * hq
  ##
  #   corrected on 05Jun2025 (needs product rule)
  ahq <- (2*a2q*t_ms_t02q + (a1q*t_ms_t0q + 2*a2q*t_ms_t0q*lvq)^2) * hq
  ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
  #    the event time part
  sat <- status*(2*a2*t_ms_t02)
  sat_p <- tapply(sat,dat$pid,sum)
  #    putting the two together
  dev2_msm <- sat_p - ahq_p

    ##   from longitudinal
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

  ##   added on 05Jun2025
  ##   Newton-Raphson to find eta such that
  ##   the first derivative is 0 (see p. 169 Rue and Held 2005)
  eta_next <- bx0 - dev1/dev2
  out <- list(eta_next=eta_next,dev1=dev1,dev2=dev2)
  return(out)
}

########################################
####  association parameters
####  (currently is using current value)
########################################
#' @export
f_association_approximated <- function(bx,bx0,params,dat,which_par='') {
    ats <- dat$transitions_to_analyse
    QD <- dat$Q
    pm0 <- params
    nms <- names(bx0)
    pm0[[which_par]][nms] <- bx0[nms]
    status <- dat$status

    dev1 <- dev2 <- rep(0,length(bx0))
    names(dev1) <- names(dev2) <- nms

    ##   hazards at quadrature points
    hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)

    ##   linear or quadratic association?
    a_od <- as.numeric(sub('a_','',which_par))
    ##   longitudinal values at event times without association parameters
    lh <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,a_od]
    ##   longitudinal values at quadrature points without association parameters
    lhq <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,a_od]

    d1 <- as.numeric(tapply(hq*lhq,dat$Q$row_id,sum))
    d2 <- as.numeric(tapply(hq*lhq^2,dat$Q$row_id,sum))
    prior_tau <- 1/100^2
    ats_to_update <- dat$model_spec$a1p
    if (which_par=='a_2') ats_to_update <- dat$model_spec$a2p
    for (k in ats_to_update) {
        ##   transition indicator
        tv <- dat$trans_indicators[,k]
        ##   first derivative
        dvll <- sum(status*tv*lh) - sum(tv*d1)
        dev1[k] <- dvll - prior_tau*bx0[k]
        ##   second derivative
        dvll <- - sum(tv*d2)
        dev2[k] <- dvll - prior_tau
    }
    ##   added on 05Jun2025
    ##   Newton-Raphson to find eta such that
    ##   the first derivative is 0 (see p. 169 Rue and Held 2005)
    eta_next <- bx0 - dev1/dev2
	out <- list(eta_next=eta_next,dev1=dev1,dev2=dev2)
    return(out)
}

########################################
####  country-level random effects in MSM
########################################
#' @export
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
  ##   added on 05Jun2025
  ##   Newton-Raphson to find eta such that
  ##   the first derivative is 0 (see p. 169 Rue and Held 2005)
  eta_next <- bx0 - dev1/dev2
  out <- list(eta_next=eta_next,dev1=dev1,dev2=dev2)
  return(out)
}

########################################
####  correlated intercept and slope
####   added on 10Jun2025
####   ********************************
####   ****  modified on 10Jul2025
####   ********************************
####    all negative diagonal entries in Q_m
####    are set to 0 (following p. 170 in RH, 
####    the line under Eq. 4.27)
####    this corresponds to "limit" the contribution from MSM data
####    when approximating that person's random effect 
####    density
########################################
#' @export
f_lmm_corr_random_effects_approximated <- function(bx,bx0,params,dat,which_par='') {
    status <- dat$status
    npds <- dat$npds
    pm0 <- params
    bx0_0 <- pm0[['u']]
    bx0_1 <- pm0[['d']]

    ##-#####################
    ##   from prior
    ##-#####################
    Q_D <- array(0,c(npds,2,2))  #  precision matrix
    Q_D[,1,1] <-  pm0[['Q_D']][1,1]
    Q_D[,2,2] <-  pm0[['Q_D']][2,2]
    Q_D[,1,2] <- Q_D[,2,1] <-  pm0[['Q_D']][1,2]
    B <- array(0,c(npds,2))      #  mean vector
    B[,1] <- pm0[['B']][1]
    B[,2] <- pm0[['B']][2]
    Qb_t <- t(c(pm0[['Q_D']][1,1]*pm0[['B']][1] + pm0[['Q_D']][1,2]*pm0[['B']][2],
                pm0[['Q_D']][2,1]*pm0[['B']][1] + pm0[['Q_D']][2,2]*pm0[['B']][2]
                ))
    Qb <- array(0,c(npds,2))
    Qb[,1] <- Qb_t[1]
    Qb[,2] <- Qb_t[2]
    m_D <- Qb

    ##-#####################
    ##   from LMM
    ##-#####################
    tau_y <- 1/pm0[['sd_y']]^2
    Q_y <- tau_y * dat$sumlong[['Q_y']]
    m_y <- tau_y * dat$sumlong[['m_y']]

    ##-#######################################
    ##   from MSM
    ##-#######################################
    ##    hazard at quadrature points
    hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
    ##    longitudinal measurements at event times without association parameters
    lv <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,1]  ##   CV
    lv2 <- compute_longitudinal_effects(pm0,dat,at_quadrature=F,include_association_par=F)[,2] ##   CV^2
    ##    longitudinal measurements at quadrature without association parameters
    lvq <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,1] ##   CV
    lvq2 <- compute_longitudinal_effects(pm0,dat,at_quadrature=T,include_association_par=F)[,2]##   CV^2
    ##   extract association parameter
	full_a1 <- full_a2 <- rep(0,dat$NTS)
    names(full_a1) <- names(full_a2) <- dat$transitions_to_analyse
    full_a1[names(pm0[['a_1']])] <- pm0[['a_1']]
    full_a2[names(pm0[['a_2']])] <- pm0[['a_2']]
    
    a1 <- full_a1[dat$jk_index]
    a1q <- full_a1[dat$Q$jk_index]
    a2 <- full_a2[dat$jk_index]
    a2q <- full_a2[dat$Q$jk_index]
    
##    a1 <- pm0[['a_1']][dat$jk_index]
##    a1q <- pm0[['a_1']][dat$Q$jk_index]
##    a2 <- pm0[['a_2']][dat$jk_index]
##    a2q <- pm0[['a_2']][dat$Q$jk_index]

    ##-#################
    ##    for intercept
    ##-#################
    t_ms_t0 <- t_ms_t02 <- t_ms_t0q <- t_ms_t02q <- 1
    ##    the quadrature part
    ahq <- (a1q*t_ms_t0q) * hq + 2*a2q*lvq*hq
    ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
    ##   the event time part
    sat <- status*(a1*t_ms_t0 + 2*a2*lv)
    sat_p <- tapply(sat,dat$pid,sum)
    ##   dev1: putting the two together
    dfdx <- sat_p - ahq_p
    ##   dev2
    ahq <- 2*a2q - (2*a2q*hq + (a1q + 2*a2q*lvq)^2 * hq)
    ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
    d2fdx2 <- -ahq_p
    ##-################
    ##    for slope
    ##-################
    ##   event times (observed transition or right censored)
    t_ms_t0 <- dat$t_ms_t0
    t_ms_t02 <- t_ms_t0^2
    ##   quadrature times
    t_ms_t0q <- dat$Q$t_ms_t0
    t_ms_t02q <- dat$Q$t_ms_t0^2
    ##    the quadrature part
    ahq <- (a1q*t_ms_t0q + 2*a2q*t_ms_t0q*lvq) * hq
    ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
    ##   the event time part
    sat <- status*(a1*t_ms_t0 + 2*a2*t_ms_t0*lv)
    sat_p <- tapply(sat,dat$pid,sum)
    ##   dev1: putting the two together
    dfdy <- sat_p - ahq_p
    ##   dev2: the quadrature part
    ahq <- 2*a2q*t_ms_t02q*hq + (a1q*t_ms_t0q + 2*a2q*t_ms_t0q*lvq)^2 * hq
    ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
    ##   dev2: the event time part
    sat <- 2*status*a2*t_ms_t02
    sat_p <- tapply(sat,dat$pid,sum)
    ##   dev2: putting the two together
    d2fdy2 <- sat_p - ahq_p

    ##-#####################
    ##   cross derivative
    ##-#####################
    ##   the quadrature part
    ahq <- 2*a2q*t_ms_t0q*hq + (a1q*t_ms_t0q + 2*a2q*t_ms_t0q*lvq)*(a1q + 2*a2q*lvq)*hq
    ##    ahq <- (a1q^2*t_ms_t0q) * hq
    ahq_p <- as.numeric(tapply(ahq,dat$Q$pid,sum))
    ##   the event time part
    sat <- 2*status*a2*t_ms_t0
    sat_p <- tapply(sat,dat$pid,sum)
    d2fdxdy <- sat_p - ahq_p

    ##-#####################
    ##   modified on 10Jul2025
    ##   p.170 in RH
    ##-#####################
##    ids <- which(d2fdx2>0 | d2fdy2>0)
##    if (length(ids)>0) {
##    	d2fdx2[ids] <- 
##    	d2fdy2[ids] <- 
##    	d2fdxdy[ids] <- 0
##    }

    ##   the Q matrix, the precision matrix, from the MSM part
    Q_m <- array(0,c(npds,2,2))
    ##    corrected on 25Jun2025
    ##    all the elements in Q should be multipled by
    ##    -1 (see joplin note on Correlated random effects
    ##    joplin//x-callback-url/openNote?id=1006cb9d2b8642b6913849739c11ed61)
    Q_m[,1,1] <- -1*d2fdx2
    Q_m[,2,2] <- -1*d2fdy2
    Q_m[,1,2] <- Q_m[,2,1] <- -1*d2fdxdy

    ##   the m vector
    m_m <- array(0,c(npds,2))
    m_m[,1] <- dfdx - bx0_0*d2fdx2 - bx0_1*d2fdxdy
    m_m[,2] <- dfdy - bx0_0*d2fdxdy - bx0_1*d2fdy2
    
    ##-#####################
    ##   modified on 10Jul2025
    ##   p.170 in RH
    ##-#####################
    V_m <- inv2x2matrix_vectorised(Q_m)
    ids <- which(V_m[,1,1]<0 | V_m[,2,2]<0)
    if (length(ids)>0) Q_m[ids,,] <- m_m[ids,] <- 0
    
    ##-########################
    ##   putting all together
    ##-########################
    Q <- V <- Q_D + Q_y + Q_m
    m <- matrix(0,nrow=2,ncol=npds)
    eta_next <- t(m)
    mm <- m_D + m_y + m_m
    m[1,] <- mm[,1]
    m[2,] <- mm[,2]
    sx <- NULL
    V <- inv2x2matrix_vectorised(Q)
    ids <- which(V[,1,1]<0 | V[,2,2]<0)
    if (length(ids)>0) {
    	for (ii in ids) {
    		##   remove contribution from the MSM data 
    		##   when constructing the proposal distribution
    		Q[ii,,] <- Q_D[ii,,] + Q_y[ii,,]
    		V[ii,,] <- MASS::ginv(Q[ii,,])
    		m[1,ii] <- m_D[ii,1] + m_y[ii,1]
    		m[2,ii] <- m_D[ii,2] + m_y[ii,2]
    	}
    }    
    eta_next <- fast_dotp(m,V)
    out <- list(eta_next=eta_next,m=m,V=V,Q=Q)
    return(out)
}

#' @export
fast_dotp <- function(m,V) {
    i1 <- m[1,]*V[,1,1] + m[2,]*V[,2,1]
    i2 <- m[1,]*V[,1,2] + m[2,]*V[,2,2]
    out <- cbind(i1,i2)
    return(out)
}

#' @export
inv2x2matrix_vectorised <- function(vm) {
    out <- vm
    out[,1,1] <- vm[,2,2]
    out[,2,2] <- vm[,1,1]
    out[,1,2] <- out[,2,1] <- -vm[,1,2]
    d <- 1/(vm[,1,1]*vm[,2,2] - vm[,1,2]^2)
    out <- d*out
    return(out)
}

#' @export
same_elements <- function(x, eps=0.001) {
  out <- T
  d <- abs((x-x[1])/x[1])
  if (any(d>eps)) out <- F
  return(out)
}

#' @export
obtain_IGMRF_Q <- function(n, s, vs=1e-10) {
#  on entry: 
#     n  = number of random effects
#     s  = (conditional) standard derivation
#     vs = a very small value assigned to the off-diagonals of the identity matrix
#          (default value = 1e-10)
    A <- matrix(vs,nrow=n,ncol=n)
    tau <- 1/s^2
    diag(A) <- tau
    E <- eigen(A)
    V <- E$vectors
    l <- E$values[2:n]
    tildeLambda <- diag(c(0,l))
    tildeQ <- V%*%tildeLambda%*%t(V)
    ##   quality control
    if (!same_elements(V[,1])) stop('not all the elements in e1 are the same')
    if (!same_elements(l)) {
        stop('not all eigenvalues are the same')
    } else {
        if (!same_elements(c(tau,l))) stop('the eigenvalues should be the same as the precision')
    }
    out <- list(tildeLambda=tildeLambda, tildeQ=tildeQ, l=l[1], V=V)
    return(out)
}

#' @export
f_w_approximated <- function(bx,bx0,params,dat,which_par='') {
  which_jk <- sub('w_','',which_par)
  QD <- dat$Q
  pm0 <- params
  nms <- names(bx0)
  pm0[['w']][nms] <- bx0[nms]
  status <- dat$status

  hq <- compute_haz(pm0,dat,at_quadrature = T,include_quadweights = T)
  hr <- as.numeric(tapply(hq,dat$Q$row_id,sum))
  
  dfdw2 <- array(0, dat$nctys)
  for (icty in 1:dat$nctys) {
    dfdw2[icty] <- sum(hr[which(dat$jkc_index==paste0(icty,'_',which_jk))])
  }
  tildeQ <- obtain_IGMRF_Q(dat$nctys,1/sqrt(pm0[['sd_w']]))$tildeQ

  A <- matrix(1,ncol=dat$nctys,nrow=1)
  E <- 0
  Q <- diag(dfdw2) + tildeQ
  ncases <- dat[['ncases_by_jkc']][paste0(1:dat$nctys,'_',which_jk)]
  m <- ncases - dfdw2 + dfdw2*pm0[['w']][paste0(1:dat$nctys,'_',which_jk)]
  m <- matrix(m,ncol=1)

  iQ <- MASS::ginv(Q)
  mu <- iQ%*%m  
  mu_star <- mu - iQ%*%t(A)%*%MASS::ginv(A%*%iQ%*%t(A))%*%(A%*%mu-E)
  s_star <- iQ - iQ%*%t(A)%*%MASS::ginv(A%*%iQ%*%t(A))%*%(A%*%iQ)

  mu_star <- mu_star[,1]
  eta_next <- mu_star
  names(eta_next) <- names(mu_star) <- nms
  out <- list(eta_next=eta_next,m=mu_star,V=s_star,Q=Q, ncases=ncases)
  return(out)
}

#' @export
logden_exchangeable_s2z <- function(x,prec) {
  n <- length(x)
  A <- matrix(1,ncol=n,nrow=1)
  E <- 0  #  sum to 0 constraint
  tau <- prec
  Q <- tau * diag(n)
  mu <- matrix(0,nrow=n,ncol=1)

  iQ <- MASS::ginv(Q)
  mu_star <- mu - iQ%*%t(A)%*%MASS::ginv(A%*%iQ%*%t(A))%*%(A%*%mu-E)
  s_star <- iQ - iQ%*%t(A)%*%MASS::ginv(A%*%iQ%*%t(A))%*%(A%*%iQ)

  ##   log density at xx
  xx <- matrix(x,nrow=n)
  ED <- eigen(s_star)
  va <- ED$values
  vas <- va
  eps <- 1e-10
  vas[which(va<eps)] <- 0
  vas[which(va>eps)] <- 1/va[which(va>eps)]
  ve <- ED$vectors
  eps <- 1e-10
  inv_s_star <- ve%*%diag(vas)%*%t(ve)
  out <- -(n-1)/2*log(2*pi) - 1/2*sum(log(va[which(va>eps)])) - 1/2*t(xx-mu_star)%*%inv_s_star%*%(xx-mu_star)
  return(as.numeric(out))
}
