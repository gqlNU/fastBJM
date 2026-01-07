# fastBJM: An efficient and scalable Bayesian algorithm for fitting joint models of longitudinal and event history data

Joint models are computationally expensive to fit due to the complex joint likelihood and the need to estimate a large number of individual-level random effects. We propose fastBJM, an efficient and scalable algorithm for fitting Bayesian joint models using Markov chain Monte Carlo.

The preprint of the paper that details the algorithm will be made avaiable on [Research Square](https://www.researchsquare.com/).


##  Installing the package

```R
devtools::install_github("gqlNU/fastBJM")
library("fastBJM")
```

The code below demonstrates the use of fastBJM to fit a simulated dataset with 5000 people (more detail to follow).

```R
#################################
##   define permissible transitions
#################################
ats <- c('12','13','15','24','25','34','35','45')
names(ats) <- ats


model_spec <- list()
model_spec$allowable_transitions <- model_spec$transitions_to_analyse <- ats
## --- name of the longitudinal variable in the dataset 
model_spec$which_longvar <- 'y'
## --- association structure:
##   aps <- c('a_1','a_2'): a quadratic function of the current value with the transition intensities
##   aps <- c('a_1')      : a linear function of the current value with the transition intensities
##   aps <- NULL          : the two submodels are not linked
model_spec$aps <- c('a_1','a_2')
model_spec$a1p <- ats
model_spec$a2p <- ats
## --- fixed effects
##   X_spec[,1]: names of the columns to be included as fixed effects
##   X_spec[,2]: types of covariates (cat=binary/categorical and cnt=continuous)
model_spec$X_spec <- matrix(c('x1','cnt',
                              'x2','cat'),byrow=T,ncol=2)
## --- which transition intensities are specified as a function of the fixed effects
##   e.g. x1=ats means x1 affects all permissible transitions and each transition has its own effect estimate
##   e.g. x1=c('12','15') means x1 only affects the 12 and 15 transitions
model_spec$beta_by_country <- list(x1=ats,
   	                               x2=ats)
## --- definition of the age intervals for the baseline intensities
##   the setting below corresponds to two age intervals, [50, 70) and [70,150]
model_spec$age_range <- c(50,150)
model_spec$age_cuts <- c(70)
##   other settings - leave as default
model_spec$nQs <- 15
model_spec$lmm_correlated_random_effects <- TRUE
model_spec$standardise_y <- T
model_spec$standardise_cnt_x <- T
model_spec$weibull_baseline <- FALSE
model_spec$joint <- !is.null(model_spec$which_longvar)
model_spec$include_fixed_effects <- !is.null(model_spec$X_spec)
model_spec$update_RE_byperson <- TRUE #  update each person independently at the MH step
model_spec$byperson <- TRUE


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
niters <- 2
sims.list <- mcmc_update(fitdata, inits, niters, model_spec, update_setting)

```

