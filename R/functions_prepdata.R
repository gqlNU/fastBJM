#' @export
prepare_fitdata <- function(dat, model_spec) {
    ##   add participant ID
    dd <- dat
    dd$msm$mergeid <- paste0('p',dd$msm[[model_spec$PID]])
    dd$long$mergeid <- paste0('p',dd$long[[model_spec$PID]])
    dat <- as.list(dd$msm)
    dat$long <- dd$long
    all_data <- dat

    ##   subset data according to transition selection by user
    out <- subset_data_by_transitions(all_data,
                                      model_spec$allowable_transitions, 
                                      model_spec$transitions_to_analyse)
    ###  checks for the real data (to work on)
    if (FALSE) {
        ##   remove people with extreme longitudinal measurements
        all_data <- subset_data_by_longitudinal(all_data,model_spec$eQs)
        ##   remove people with extremely large/small longitudinal predictions
        if (is.null(simdata)) all_data <- subset_data_by_sensible_longitudinal(all_data,model_spec$eQs)
        ##   remove people whose last recorded age is far away from their last longitudinal measurement
        all_data <- subset_data_by_last_age(all_data, model_spec$max_diff)
        ##   remove people whose age of entry is far away from their first longitudinal measurement
        all_data <- subset_data_by_first_age(all_data, model_spec$max_diff)
        out <- all_data
    }
    ##   prepare fixed effects if included
    if (!is.null(model_spec$X_spec)) out <- prepare_X(out, model_spec)
    ##   prepare age-related variables
    out <- prepare_age(out, model_spec)
    ##   prepare age-related variables
    out <- prepare_indices(out)
    ##   get GMRF settings
    out$gmrf_settings <- gmrf_configs()
    ##   for Gauss Legendre qradrature
    out$nQs <- model_spec$nQs
    out$Q <- prepare_quadrature(out)
    ##   unique subject IDs (use this throughout!!)
    out$subject <- unique(out$mergeid)
    ##   structure and summarise longitudinal measurements
    which_longvar <- model_spec$which_longvar
    if (!is.null(which_longvar)) {
        out$long <- format_longitudinal(out, model_spec)
        out$sumlong <- summarise_longitudinal_measures(out)
    }
    out$model_spec <- model_spec
    print(out$ncases_by_jkg)
    return(out)
}


########################################
####  subset data by user-selected
####  transitions
########################################
#' @export
subset_data_by_transitions <- function(dat, allowable_transitions, transitions_to_analyse) {
    lg <- klg <- dat$long
    if (any(names(dat)=='msm')) {
        ddt <- dat$msm
    } else {
        ddt <- as.data.frame(dat[-c(which(names(dat)=='long'))])
    }
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
    ##   get the longitudinal measurements of the people we keep
    pds <- unique(out$mergeid)
    ids <- NULL
    for (pd in pds) {
        if (length(which(lg$mergeid==pd))>0) {
            ids <- c(ids, which(lg$mergeid==pd))
        }
    }
    klg <- lg[ids,]
    out$long <- klg
    out$allowable_transitions <- allowable_transitions
    out$transitions_to_analyse <- transitions_to_analyse
    out$NTS <- length(out$transitions_to_analyse)
    return(out)
}

#######################################################
###  function to get fixed effects
#######################################################
#' @export
prepare_X <- function(dat, model_spec) {
  out <- dat
  X_spec <- model_spec$X_spec
  out$X_spec <- data.frame(xn=X_spec[,1],type=X_spec[,2])
  out$X <- gather_X(out, model_spec$standardise_cnt_x)
  out$nXs <- ncol(out$X)
  out$XVars <- colnames(out$X)
  return(out)
}


#######################################################
###  function to check the validity of the fixed
###  effect specification
#######################################################
#' @export
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
#' @export
gather_X <- function(dat,standardise_cnt_x=T) {
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
            if (length(unique(x))>10) print(paste0('Warning: there are more than 10 levels in the categorical covariate ',this_x))
            xc <- fastDummies::dummy_cols(x,remove_first_dummy=TRUE)
            xc <- as.matrix(xc[which(names(xc)!='.data')])
            this_xc <- paste0(this_x,1:ncol(xc))
            out <- cbind(out,xc)
            out_colnames <- c(out_colnames,this_xc)
        }
        if (tp=='cnt') {
            #  a continuous one
            xs <- x
            if (standardise_cnt_x) xs <- scale(x,center=TRUE,scale=TRUE)
            out <- cbind(out,xs)
            out_colnames <- c(out_colnames,this_x)
        }
    }
    out <- as.matrix(out)
    colnames(out) <- out_colnames
  }
  return(out)
}


#######################################################
###  function to define age intervals
#######################################################
#' @export
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
#' @export
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
#' @export
prepare_age <- function(dat, model_spec) {
  out <- dat
  age_range <- model_spec$age_range
  age_cuts <- model_spec$age_cuts
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
#' @export
prepare_indices <- function(dat) {
    out <- dat
    ats <- out$transitions_to_analyse
    out$unique_mergeid <- unique(out$mergeid)
    out$npds <- length(out$unique_mergeid)
    out$unique_mergeid_id <- 1:out$npds
    names(out$unique_mergeid_id) <- out$unique_mergeid
    out$pid <- out$unique_mergeid_id[as.character(out$mergeid)]

    ##   country-related
    if (any(names(dat)=='country')) {
        out$cty_id <- as.numeric(as.factor(out$country))
    } else {
        out$country <- 1
        out$cty_id <- as.numeric(as.factor(out$country))
    }
    names(out$cty_id) <- out$country
    out$nctys <- length(unique(out$cty_id))

    ##   index of transition
    out$jk_index <- paste0(out$from,out$to)
    ##   index of transition-country
    out$jkc_index <- paste0('cty',out$cty_id,'_',out$jk_index)
    ##   index of transition-age
    out$jkg_index <- paste0(out$jk_index,'g',out[['age_gp_calc']])
    ##   index of transition-age_country
    out$jkgc_index <- paste0('cty',out$cty_id,'_',out$jk_index,'g',out[['age_gp_calc']])

    out$all_jkg_index <- c(sapply(ats,function(k){paste0(k,'g',1:out$nGs)}))
    out$ncases_by_jkg <- rep(0,length(out$all_jkg_index))
    names(out$ncases_by_jkg) <- out$all_jkg_index
    ts <- tapply(out$status,out$jkg_index,sum)
    out$ncases_by_jkg[names(ts)] <- ts[names(ts)]

    out$trans_indicators <- array(0,c(length(out$from),out$NTS))
    colnames(out$trans_indicators) <- out$transitions_to_analyse
    for (k in colnames(out$trans_indicators))
        out$trans_indicators[which(out$jk_index==k),k] <- 1

    return(out)
}

########################################
####  GMRF related configs
########################################
#' @export
gmrf_configs <- function() {
    out <- list()
    #    more search points improve acceptance at the cost of a slightly longer run time
    out$search_range <- c(-4,4)
    out$search_points <- 1000
    #  auto_mode_finder
    out$eps <- 5        #  accruacy for the auto mode finder abs((new-old)/new)*100 < eps
    out$max_iters <- 50 #  max iterations allowed for auto mode finder
    return(out)
}


########################################
####  for Gauss Legendre quadrature
########################################
#' @export
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
  Q$mergeid <- rep(dat$mergeid,each=nQs)
  Q$row_id <- rep(1:length(dat$pid),each=nQs)
  Q$age_entry <- rep(dat$age_entry,each=nQs)
  Q$t_ms_t0 <- Q$age_calc - Q$age_entry

  return(Q)
}



#######################################################
###  function to format the longitudinal measurements
#######################################################
#' @export
format_longitudinal <- function(dat, model_spec) {
    which_longvar <- model_spec$which_longvar
    standardise_y <- model_spec$standardise_y
    dlong <- dat$long
    dlong$pid <- dat$unique_mergeid_id[as.character(dlong$mergeid)]
    ##   time of longitudinal measurement since entering SHARE
    tmp <- dat$age_entry
    names(tmp) <- dat$mergeid
    if (length(grep('age_entry',names(dlong)))==0) {
        dlong$age_entry <- tmp[as.character(dlong$mergeid)]
    }
    dlong$t_ms_t0 <- dlong$age - dlong$age_entry
    ##   extract the longitudinal variable
    y <- dlong[[which_longvar]]
    mm <- mean(y)
    ss <- sd(y)
    if (standardise_y) y <- (y - mm)/ss
    dlong$y <- y
    dlong <- as.list(dlong)
    dlong$nts <- table(dlong$pid)
    dlong$mean_y <- mm
    dlong$sd_y <- ss
    return(dlong)
}


#######################################################
###  convert output from tapply to a data frame
#######################################################
#' @export
to_dataframe <- function(x,nm) {
    out <- data.frame(mergeid=names(x))
    out[[nm]] <- x
    return(out)
}

#######################################################
###  summarise longitudinal data
#######################################################
#' @export
summarise_longitudinal_measures <- function(dat) {
    lg <- dat$long
    n <- to_dataframe(tapply(lg$y,lg$mergeid,length),nm='n')
    h <- lg$age - lg$age_entry  #  years since entry
    ##   various summaries
    sum_h <- to_dataframe(tapply(h,lg$mergeid,sum),nm='sum_h')
    sum_h2 <- to_dataframe(tapply(h^2,lg$mergeid,sum),nm='sum_h2')
    sum_y <- to_dataframe(tapply(lg$y,lg$mergeid,sum),nm='sum_y')
    sum_hy <- to_dataframe(tapply(h*lg$y,lg$mergeid,sum),nm='sum_hy')
    out <- merge(n,sum_h,by='mergeid',sort=F)
    out <- merge(out,sum_h2,by='mergeid',sort=F)
    out <- merge(out,sum_y,by='mergeid',sort=F)
    out <- merge(out,sum_hy,by='mergeid',sort=F)
    ##   reorder out according to dat$subject
    ids <- sapply(dat$subject,function(x){which(out$mergeid==x)})
    out <- out[ids,]
    npds <- length(n$mergeid)
    Q_y <- array(0,c(npds,2,2))
    Q_y[,1,1] <- out$n
    Q_y[,2,2] <- out$sum_h2
    Q_y[,1,2] <- Q_y[,2,1] <- out$sum_h
    m_y <- array(0,c(npds,2))
    m_y[,1] <- out$sum_y
    m_y[,2] <- out$sum_hy
    dimnames(Q_y)[[1]] <- dat$subject
    dimnames(m_y)[[1]] <- dat$subject
    out <- as.list(out)
    out$Q_y <- Q_y
    out$m_y <- m_y
    return(out)
}