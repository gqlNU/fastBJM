prepare_fitdata <- function(icty, gender, data_type,
                            allowable_transitions, transitions_to_analyse,
                            X_spec,
                            age_range, age_cuts,
                            nQs,
                            which_longvar) {
  #  define which countries to be included
  all_ctys <- country_lookup()
  ctys <- all_ctys$country
  if (!is.null(icty)) ctys <- all_ctys$country[icty]  #  how many countries to be included in this fit
  #if (!is.null(icty)) print(all_ctys[icty,])

  #  read all transition data from the selected countries
  all_data <- read_final_data(ctys, gender, type=data_type)
  #  subset data according to transition selection by user
  out <- subset_data_by_transitions(all_data, allowable_transitions, transitions_to_analyse)
  #  prepare fixed effects if included
  out <- prepare_X(out, X_spec)
  #  prepare age-related variables
  out <- prepare_age(out, age_range, age_cuts)
  #  prepare age-related variables
  out <- prepare_indices(out)
  #  get GMRF settings
  out$gmrf_settings <- gmrf_configs()
  #  for Gauss Legendre qradrature
  out$nQs <- nQs
  out$Q <- prepare_quadrature(out)
  #  structure longitudinal measurements
  out$long <- format_longitudinal(out, which_longvar)
  return(out)
}


########################################
####  subset data by user-selected
####  transitions
########################################
subset_data_by_transitions <- function(dat, allowable_transitions, transitions_to_analyse) {
  lg <- dat$long
  ddt <- as.data.frame(dat[-c(which(names(dat)=='long'))])
  trjs <- paste0(ddt[['from']],ddt[['to']])
  jis <- rep(0,length(trjs))
  for (k in transitions_to_analyse) {
    ii <- which(trjs==k)
    if (length(ii)>0) jis[ii] <- 1
  }
  ddt <- ddt[which(jis==1),]
  ddt <- as.list(ddt)
  trjs <- paste0(ddt[['from']],ddt[['to']])
  print(table(trjs,ddt$status))
  out <- ddt
  out$long <- lg
  out$allowable_transitions <- allowable_transitions
  out$transitions_to_analyse <- transitions_to_analyse
  out$NTS <- length(out$transitions_to_analyse)
  return(out)
}

########################################
####  GMRF related configs
########################################
gmrf_configs <- function() {
    out <- list()
    #    more search points improve acceptance at the cost of a slightly longer run time
    out$search_range <- c(-4,4)
    out$search_points <- 1000
    #  auto_mode_finder
    out$eps <- 1  		#  accruacy for the auto mode finder abs((new-old)/new)*100 < eps
    out$max_iters <- 20	#  max iterations allowed for auto mode finder
    return(out)
}


########################################
####  for Gauss Legendre quadrature
########################################
prepare_quadrature <- function(dat) {
  nQs <- dat$nQs
  glq <- statmod::gauss.quad(nQs, kind = "legendre")
  #   constructing Q
  Q <- list()
  nrs <- length(dat[['start']])
  sojourn <- rep(dat[['stop']] - dat[['start']],each=nQs) # the duration at state j within the jk transition
  Q$w <- sojourn/2*rep(glq$weights,nrs)

  tts <- rep((glq$nodes+1)/2,nrs)  # change the [-1,1] support to [0,1]
  stts <- sojourn * tts            # change the [0,1] support to [0,sojourn] for each transition
  start_age <- rep(dat[['start']],each=nQs)  #  age at the start of each transition
  Q$age_calc <- start_age + stts
  Q$age_gp_calc <- which_age_interval_cut(dat,Q$age_calc)

  #  various indices
  q_from <- rep(dat$from,each=nQs)
  q_to <- rep(dat$to,each=nQs)
  Q$jk_index <- paste0(q_from,q_to)
  Q$jkg_index <- paste0(Q$jk_index,'g',Q$age_gp_calc)
  Q$jkc_index <- rep(dat$jkc_index,each=nQs)
  Q$jkgc_index <- rep(dat$jkgc_index,each=nQs)

  #  country index
  Q$cty_id <- rep(dat$cty_id,each=nQs)
  Q$pid <- rep(dat$pid,each=nQs)
  Q$row_id <- rep(1:length(dat$pid),each=nQs)
  Q$age_entry <- rep(dat$age_entry,each=nQs)
  Q$t_ms_t0 <- Q$age_calc - Q$age_entry

  return(Q)
}

#######################################################
###  function to define age intervals
#######################################################
get_age_intervals <- function(age_range,age_cuts) {
    if (!is.null(age_cuts)) {
    al <- c(age_range[1],age_cuts)
    au <- c(age_cuts,age_range[2])
    } else {
    al <- age_range[1]
    au <- age_range[2]
    }
    age_intervals <- cbind(al,au)
    colnames(age_intervals) <- c('lower','upper')
    rownames(age_intervals) <- paste0(age_intervals[,1],'-',age_intervals[,2])
    return(age_intervals)
}

#######################################################
###  allocating ages to age groups
#######################################################
which_age_interval_cut <- function(dat, ages) {
  bks <- unique(c(dat$age_intervals))
  ncuts <- length(bks)-1
  cts_nms <- cut(ages,bks,right=F)
  out <- as.numeric(cut(ages,bks,labels=1:ncuts,right=F))
  names(out) <- cts_nms
  return(out)
}

#######################################################
###  gather various age-related variables
#######################################################
prepare_age <- function(dat, age_range, age_cuts) {
  out <- dat
  #  age-related
  age_intervals <- get_age_intervals(age_range, age_cuts)
  out$age_intervals <- age_intervals
  out$age_range <- age_range
  out$age_cuts <- age_cuts
  out$nGs <- nrow(out$age_intervals)
  out$age_calc <- out$stop  #  age at the end of the stay at a state
  out$t_ms_t0 <- out$stop - out$age_entry  #  event time since entry
  out$age_gp_calc <- which_age_interval_cut(out, out[['age_calc']])
  return(out)
}

#######################################################
###  gather various indices
#######################################################
prepare_indices <- function(dat) {
  out <- dat
  ats <- out$transitions_to_analyse
  out$unique_mergeid <- unique(out$mergeid)
  out$npds <- length(out$unique_mergeid)
  out$unique_mergeid_id <- 1:out$npds
  names(out$unique_mergeid_id) <- out$unique_mergeid
  out$pid <- out$unique_mergeid_id[out$mergeid]

#  country-related
  out$cty_id <- as.numeric(as.factor(out$country))
  names(out$cty_id) <- out$country
  out$nctys <- length(unique(out$cty_id))

  out$jk_index <- paste0(out$from,out$to)
  out$jkc_index <- paste0('cty',out$cty_id,'_',out$jk_index)
  out$jkg_index <- paste0(out$jk_index,'g',out[['age_gp_calc']])
  out$jkgc_index <- paste0('cty',out$cty_id,'_',out$jk_index,'g',out[['age_gp_calc']])

  out$all_jkg_index <- c(sapply(ats,function(k){paste0(k,'g',1:out$nGs)}))
  out$ncases_by_jkg <- rep(0,length(out$all_jkg_index))
  names(out$ncases_by_jkg) <- out$all_jkg_index
  ts <- tapply(out$status,out$jkg_index,sum)
  out$ncases_by_jkg[names(ts)] <- ts[names(ts)]

  out$trans_indicators <- array(0,c(length(out$from),out$NTS))
  colnames(out$trans_indicators) <- out$transitions_to_analyse
  for (k in colnames(out$trans_indicators)) out$trans_indicators[which(out$jk_index==k),k] <- 1
  
  return(out)
}

#######################################################
###  function to calculate the baseline hazards
#######################################################
compute_baseline <- function(params, dat, at_quadrature) {
  #  age- and transition-specific baseline hazards
  dd <- dat
  if (at_quadrature) dd <- dat$Q
  h0 <- params[['l']][dd[['jkg_index']]]
  out <- h0
  return(out)
}

#######################################################
###  function to get fixed effects
#######################################################
prepare_X <- function(dat,X_spec) {
  out <- dat
  if (!is.null(X_spec)) {
    out$X_spec <- data.frame(xn=X_spec[,1],type=X_spec[,2])
    out$X <- gather_X(out)
    out$nXs <- ncol(out$X)
    out$XVars <- colnames(out$X)
  }
  return(out)
}

#######################################################
###  function to check the validity of the fixed
###  effect specification
#######################################################
check_X_spec <- function(dat) {
    out <- NULL
    Xs <- dat$X_spec
    for (ix in 1:nrow(Xs)) {
        this_x <- Xs[['xn']][ix]
        tp <- Xs[['type']][ix]
        err <- ''
        if (all(names(dat)!=this_x)) 
          err <- paste0(this_x,' is not in the dataset')
        if (err!='') out <- c(out,err)
        err <- ''
        if (all(c('cnt','cat')!=tp)) 
          err <- paste0('type of ',this_x,' must be cat or cnt')
        if (err!='') out <- c(out,err)
    }
    if (!is.null(out)) {
        print('Error in X_spec')
        print(out)
        stop('------')
    }
    return(out)
}

#######################################################
###  function to construct the X matrix
#######################################################
gather_X <- function(dat) {
  Xs <- dat$X_spec
  out <- out_colnames <- NULL
  if (!is.null(Xs)) {
    #  check if all is OK in X_spec
    err <- check_X_spec(dat)
    #  go through each fixed effect
    for (ix in 1:nrow(Xs)) {
        this_x <- Xs[['xn']][ix]
        tp <- Xs[['type']][ix]
        x <- dat[[this_x]]
        if (tp=='cat') {
            #  a categorical one
            if (length(unique(x)>10)) print(paste0('Warning: there are more than 10 levels in the categorical covariate ',this_x))
            xc <- fastDummies::dummy_cols(x,remove_first_dummy=TRUE)
            xc <- as.matrix(xc[which(names(xc)!='.data')])
            this_xc <- paste0(this_x,1:ncol(xc))
            out <- cbind(out,xc)
            out_colnames <- c(out_colnames,this_xc)
        }
        if (tp=='cnt') {
            #  a continuous one
            out <- cbind(out,x)
            out_colnames <- c(out_colnames,this_x)
        }
    }
    out <- as.matrix(out)
    colnames(out) <- out_colnames
  }
  return(out)
}


#######################################################
###  function to format the longitudinal measurements
#######################################################
format_longitudinal <- function(dat, which_longvar='bmi') {
  dlong <- dat$long
  dlong$pid <- dat$unique_mergeid_id[dlong$mergeid]
  #  time of longitudinal measurement since entering SHARE
  tmp <- dat$age_entry
  names(tmp) <- dat$mergeid
  dlong$age_entry <- tmp[dlong$mergeid]
  dlong$t_ms_t0 <- dlong$age - dlong$age_entry
  #  standardise longitudinal measurements (y-mean(y))/sd(y)
  y <- dlong[[which_longvar]]
  dlong$y <- (y - mean(y))/sd(y)
  dlong <- as.list(dlong)
  dlong$nts <- table(dlong$pid)
  return(dlong)
}
