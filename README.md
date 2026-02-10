# fastBJM: An efficient and scalable Bayesian algorithm for fitting joint models of longitudinal and event history data

Joint models are computationally expensive to fit due to the complex joint likelihood and the need to estimate a large number of individual-level random effects. `fastBJM` is an efficient and scalable algorithm for fitting Bayesian joint models using Markov chain Monte Carlo.

The preprint of the paper that details the algorithm will be made avaiable on [Research Square](https://www.researchsquare.com/).


##  Installing the package

```R
devtools::install_github("gqlNU/fastBJM")
```

## Example 1: Fitting a joint model

This example illustrates the use of fastBJM to analyse a simulated dataset with the longitudinal and event history data from 5000 people. The multistate submodel is given by the state diagram below. There are five health states (1 - 5) and eight permissible transitions $\Delta=\{(12), (13), (15), (24), (35), (34), (35), (45)\}$. For each person, the first longitudinal measurement is taken at the start of the follow-up and subsequent measurements are taken every 2.5 years until the person exits the study due to either (non-informative) drop-out or death (which is State 5 in the multistate model).

<img width="601" height="467" alt="Screenshot 2026-01-08 at 15 17 29" src="https://github.com/user-attachments/assets/48cfe589-99f5-4913-9102-b1ef2639257a" />

The example code below fits a joint model in which a linear growth curve model with random individual-specific intercepts and slopes is used for the longitudinal submodel. For the multistate submodel, the transition intensities are specified as functions of the effects from two (fixed-effect) covariates, $x_1$ (continuous) and $x_2$ (binary) and a linear association structure with the current value of the underlying longitudinal process. The intensity of a generic permissible transition from State $j$ to State $k$ for person $i$ takes the following structure

$$
\lambda^{(jk)}_i(t)=\lambda^{(jk)}_0(t)\cdot e^{\boldsymbol{X}_i\boldsymbol{\beta}^{(jk)} + a_1(u_i+d_i t)}
$$

where $\boldsymbol{X}_i$ denotes the covariates and $u_i$ and $d_i$ are the intercept and slope. The code below specifies a piecewise constant model on the baseline intensity $\lambda^{(jk)}_0(t)$ with two age-intervals 50-70 and 70-150.

<details>
<summary> R code for fitting a joint model </summary>

```R
library("fastBJM")

#################################
##   user inputs
#################################
##   define permissible transitions
ats <- c('12','13','15','24','25','34','35','45')
names(ats) <- ats
model_spec <- list()
##   person ID that links the event history data to the longitudinal data
model_spec$PID <- 'mergeid'
##   the set of transitions to be included in the analysis
model_spec$allowable_transitions <- model_spec$transitions_to_analyse <- ats
## --- name of the longitudinal variable in the dataset
model_spec$which_longvar <- 'y'
## --- association structure:
##   aps <- c('a_1','a_2'): a quadratic function of the current value with the transition intensities
##   aps <- c('a_1')      : a linear function of the current value with the transition intensities
##   aps <- NULL          : for a separate analysis of the two submodels
model_spec$aps <- c('a_1','a_2')
model_spec$a1p <- ats
model_spec$a2p <- ats
## --- fixed effects on the transition intensities
##   X_spec[,1]: names of the columns to be included as fixed effects
##   X_spec[,2]: types of covariates (cat=binary/categorical and cnt=continuous)
##   set model_spec$X_spec <- NULL if no fixed effect included
model_spec$X_spec <- matrix(c('x1','cnt',
                              'x2','cat'),byrow=T,ncol=2)
## --- which transition intensities are specified as a function of the fixed effects
##   e.g. x1=ats means x1 affects all permissible transitions and each transition has its own effect estimate
##   e.g. x1=c('12','15') means x1 only affects the 12 and 15 transitions
##   beta_by_transitions will be automatically set to NULL if model_spec$X_spec is NULL (i.e. no fixed effect)
model_spec$beta_by_transitions <- list(x1=ats,
   	                                   x2=ats)
## --- definition of the age intervals for the baseline intensities
##   the setting below corresponds to two age intervals, [50, 70) and [70,150]
model_spec$age_range <- c(50,150)
model_spec$age_cuts <- c(70)

## --- include transition specific IID random effects?
model_spec$include_msm_random <- FALSE
#################################
##   end user input
#################################

##   specify other specifications (using defaults)
model_spec <- tidy_model_spec(model_spec)

##   get data in
dfile <- system.file("extdata", "simdata1.rds", package = "fastBJM")
dd <- readRDS(dfile)
fitdata <- prepare_fitdata(dd,model_spec)

##   define updating scheme based on the model specification
update_setting <- gather_update_setting(fitdata, model_spec)

##   specify model parameters and their starting values
params <- get_parameters(fitdata, model_spec)
inits <- initialise_parameters(fitdata, params, model_spec, update_setting)

##   carry out MCMC update
niters <- 10
sims.list <- mcmc_update(fitdata, inits, niters, model_spec, update_setting)

```
</details>



## Example 2: Fitting a multistate only model

fastBJM can be used to fit a multistate-only model with a piecewise constant model on the transition intensities. This is done by not specifying a longitudinal variable in the model specification, namely, by setting `model_spec$which_longvar <- NULL`. Fixed effects can still be included on all (or some) of the transition intensities.

To illustrate, the code below fits the six-states model for leukemia patients after bone marrow transplantation as described in [de Wreede et al. (2011)](https://www.jstatsoft.org/article/view/v038i07) using both `mstate` and `fastBJM`. For the `fastBJM` fit, the baseline transitions are time invariant and this is equivalent to having an exponential distribution in the `mstate` setting.

<details>
<summary> R code to fit the six-state model using mstate and fastBJM</summary>

```R
library('mstate')
library('fastBJM')

################################################################################
##   the code for mstate fit
##   from https://www.r-bloggers.com/2018/12/simulating-multi-state-models-with-r/
################################################################################
tmat <- mstate::transMat(x = list(c(2, 3, 5, 6), 
                         c(4, 5, 6), 
                         c(4, 5, 6), 
                         c(5, 6),
                         c(),
                         c()),
                       names = c("Tx", "Rec", "AE", "Rec+AE", 
                                 "Rel", "Death"))
data("ebmt4")
msebmt <- msprep(data = ebmt4, trans = tmat, 
                 time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"), 
                 keep = c("match", "proph", "year", "agecl"))
msebmt[msebmt$id == 2, ]

library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE)
fits_exp <- vector(mode = "list", length = n_trans)
msebmt$years <- msebmt$time/365.25
for (i in 1:n_trans){
  fits_exp[[i]] <- flexsurvreg(Surv(years, status) ~ match ,
                       data = subset(msebmt, trans == i),
                       dist = "exp")
}
jk <- unique(cbind(paste0(msebmt$from,msebmt$to),msebmt$trans))


################################################################################
##   the code for fastBJM fit
################################################################################
##   user input
####################
##   define permissible transitions
ats <- c('12','13','15','16','24','25','26','34','35','36','45','46')
names(ats) <- ats
model_spec <- list()
model_spec$allowable_transitions <- ats
model_spec$transitions_to_analyse <- ats
## --- name of the longitudinal variable in the dataset 
model_spec$which_longvar <- NULL
## --- association structure:
##   aps <- c('a_1','a_2'): a quadratic function of the current value with the transition intensities
##   aps <- c('a_1')      : a linear function of the current value with the transition intensities
##   aps <- NULL          : the two submodels are not linked
model_spec$aps <- NULL
model_spec$a1p <- ats
model_spec$a2p <- ats
## --- fixed effects
##   X_spec[,1]: names of the columns to be included as fixed effects
##   X_spec[,2]: types of covariates (cat=binary/categorical and cnt=continuous)
##   set model_spec$X_spec <- NULL if no fixed effect included
model_spec$X_spec <- matrix(c('match','cat'),byrow=T,ncol=2)
## --- which transition intensities are specified as a function of the fixed effects
##   e.g. x1=ats means x1 affects all permissible transitions and each transition has its own effect estimate
##   e.g. x1=c('12','15') means x1 only affects the 12 and 15 transitions
##   beta_by_transitions will be automatically set to NULL if model_spec$X_spec is NULL (i.e. no fixed effect)
model_spec$beta_by_transitions <- list(match=model_spec$transitions_to_analyse)
## --- definition of the age intervals for the baseline intensities
##   the setting below corresponds to two age intervals, [50, 70) and [70,150]
model_spec$age_range <- c(0,1500)
model_spec$age_cuts <- NULL
##   person ID that links the event history data to the longitudinal data
model_spec$PID <- 'id'

## --- include transition specific IID random effects?
model_spec$include_msm_random <- FALSE
####################
##   end user input
####################

##  == specify other specifications (using defaults)
model_spec <- tidy_model_spec(model_spec)

##  == get data in
msm <- msebmt[c('id','from','to','Tstart','Tstop','status','match')]
msm$match <- as.numeric(msm$match)-1
msm$age_entry <- 0
names(msm)[which(names(msm)=='Tstart')] <- 'start'
names(msm)[which(names(msm)=='Tstop')] <- 'stop'
msm$start <- msm$start/365.25
msm$stop <- msm$stop/365.25
##  == longitudinal measurements are not modelled but still need to have something in its place
long <- data.frame(id=unique(msm$id))
long$y <- rnorm(length(long))
dd <- list(msm=msm,long=long)
fitdata <- prepare_fitdata(dd,model_spec)

##  == define updating scheme based on the model specification
update_setting <- gather_update_setting(fitdata, model_spec)

##  == specify model parameters and their starting values
params <- get_parameters(fitdata, model_spec)
inits <- initialise_parameters(fitdata, params, model_spec, update_setting)

##  == carry out MCMC update
niters <- 100
res <- mcmc_update(fitdata, inits, niters, model_spec, update_setting)
sims.list <- res$sims.list


#################################
##   compare two fits
#################################
posterior_summary <- function(x) c(mean(x),quantile(x,c(0.025,0.975)))
par(mfcol=c(4,ceiling(length(model_spec$transitions_to_analyse)/2)))
ic <- 0
for (which_jk in model_spec$transitions_to_analyse) {
    ic <- ic + 1
    which_fit <- as.numeric(jk[which(jk[,1]==which_jk),2])
    res <- fits_exp[[which_fit]]$res
    plot(sims.list$l[,ic],main=paste0('l ',which_jk))
    abline(h=res['rate',c(1,2,3)],col='orange',lty=c(1,2,2),lwd=3)
    b <- NULL
    if (model_spec$include_fixed_effects) {
        plot(sims.list$beta[,ic],main=paste0('beta ',which_jk))
        abline(h=res['matchgender mismatch',c(1,2,3)],col='orange',lty=c(1,2,2),lwd=3)
        b <- round(posterior_summary(sims.list$beta[,ic]),digits=3)
    }
    print(' ==================================================================== ')
    print(paste0('--- transition (',which_jk,'): first block = fastBJM; second block = mstate'))
    print(rbind(round(posterior_summary(sims.list$l[,ic]),digits=3),b)) 
    print(round(res[,c(1,2,3)],digits=3))
}
```
</details>


## On data structure:
The data required from the user is a list that consists of two dataframes, `msm` and `long`, containing the event history data and the longitudinal measurements, respectively. The package contains a simulated dataset that illustrates the contents of the data list:

```R
dfile <- system.file("extdata", "simdata1.rds", package = "fastBJM")
dd <- readRDS(dfile)
```

To understand the data structure, we take person 4054 as an example.
```r
######  the event history data
> dd$msm[which(dd$msm$mergeid==4054),]
         start     stop from to status age_entry    x1 x2 mergeid
13391 64.00000 65.09222    1  2      1        64 -0.65  1    4054
23391 64.00000 65.09222    1  3      0        64 -0.65  1    4054
33391 64.00000 65.09222    1  5      0        64 -0.65  1    4054
41422 65.09222 73.00000    2  4      0        64 -0.65  1    4054
51422 65.09222 73.00000    2  5      0        64 -0.65  1    4054
```
This person entered the study at age 64, was diagnosed with hypertension at age 65.09 then exited at age 73 without experiencing any other health outcomes that are of interest in the analysis. Therefore, the first three rows of the `msm` data tell us that out of the three possible transitions from state 1, namely (12), (13) and (15), only the (12) transition was observed (thus with `status` being 1 on that row) while the other two with status set at 0. The subsequent two rows in `msm` tell us that this person did not experience the (24) or the (25) transition and exited at age 73. The two columns, `x1` and `x2`, are the two time-invariant covariates that we have for the individual. This is the same structure as described on p. 6 in [de Wreede et al. (2011)](https://www.jstatsoft.org/article/view/v038i07).

The longitudinal measurements for this person are stored under the column `y` in `dd$long`. The variable `age` shows the age at which each measurement was obtained:
```r
######  the longitudinal measurements
> dd$long[which(dd$long$mergeid==4054),]
      mergeid age          y age_entry
13297    4054  64 -0.8176417        64
23297    4054  66 -0.2077799        64
33176    4054  68 -0.3878730        64
42848    4054  70 -1.3820206        64
52005    4054  72 -1.0251801        64
```


## Example 3: A joint model with exchangeable random effects on each MSM transition <font color='orange'>(for illustration only)</font>
<details>
<summary> R code to fit a five-state joint model with exchangeable random effects on each transition</summary>

```R
library('fastBJM')

#################################
##   user inputs
#################################
##   define permissible transitions
ats <- c('12','13','15','24','25','34','35','45')
names(ats) <- ats
model_spec <- list()
##   person ID that links the event history data to the longitudinal data
model_spec$PID <- 'mergeid'
##   the set of transitions to be included in the analysis
model_spec$allowable_transitions <- model_spec$transitions_to_analyse <- ats
## --- name of the longitudinal variable in the dataset
model_spec$which_longvar <- 'y'
## --- association structure:
##   aps <- c('a_1','a_2'): a quadratic function of the current value with the transition intensities
##   aps <- c('a_1')      : a linear function of the current value with the transition intensities
##   aps <- NULL          : for a separate analysis of the two submodels
model_spec$aps <- c('a_1','a_2')
model_spec$a1p <- ats
model_spec$a2p <- ats
## --- fixed effects on the transition intensities
##   X_spec[,1]: names of the columns to be included as fixed effects
##   X_spec[,2]: types of covariates (cat=binary/categorical and cnt=continuous)
##   set model_spec$X_spec <- NULL if no fixed effect included
model_spec$X_spec <- matrix(c('x1','cnt',
                              'x2','cat'),byrow=T,ncol=2)
## --- which transition intensities are specified as a function of the fixed effects
##   e.g. x1=ats means x1 affects all permissible transitions and each transition has its own effect estimate
##   e.g. x1=c('12','15') means x1 only affects the 12 and 15 transitions
##   beta_by_transitions will be automatically set to NULL if model_spec$X_spec is NULL (i.e. no fixed effect)
model_spec$beta_by_transitions <- list(x1=ats,
                                       x2=ats)
## --- definition of the age intervals for the baseline intensities
##   the setting below corresponds to two age intervals, [50, 70) and [70,150]
model_spec$age_range <- c(50,150)
model_spec$age_cuts <- NULL

## --- include transition specific IID random effects?
model_spec$include_msm_random <- TRUE
jump_sd <- list()
jump_sd[['sd_w']] <- 0.2
model_spec$jump_sd <- jump_sd
model_spec$transitions_with_random <- ats
#################################
##   end user input
#################################

##   specify other specifications (using defaults)
model_spec <- tidy_model_spec(model_spec)

##   get data in
dfile <- system.file("extdata", "simdata1.rds", package = "fastBJM")
dd <- readRDS(dfile)

##   add a country indicator in (for illustration only)
nps <- length(unique(dd$msm$mergeid))
set.seed(10)
cty <- sample(1:6,nps,replace=T)
names(cty) <- unique(dd$msm$mergeid)
dd$msm$country <- dd$msm$ctyn <- cty[as.character(dd$msm$mergeid)]

fitdata <- prepare_fitdata(dd,model_spec)

#   define updating scheme based on the model specification
update_setting <- gather_update_setting(fitdata, model_spec)

##   specify model parameters and their starting values
params <- get_parameters(fitdata, model_spec)
inits <- initialise_parameters(fitdata, params, model_spec, update_setting)


##   carry out MCMC update (for illustration only)
niters <- 2
res <- mcmc_update(fitdata, inits, niters, model_spec, update_setting)
sims.list <- res$sims.list
current_pars <- res$current_pars
```
</details>
