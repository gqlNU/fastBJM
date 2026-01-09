# fastBJM: An efficient and scalable Bayesian algorithm for fitting joint models of longitudinal and event history data

Joint models are computationally expensive to fit due to the complex joint likelihood and the need to estimate a large number of individual-level random effects. We propose fastBJM, an efficient and scalable algorithm for fitting Bayesian joint models using Markov chain Monte Carlo.

The preprint of the paper that details the algorithm will be made avaiable on [Research Square](https://www.researchsquare.com/).


##  Installing the package

```R
devtools::install_github("gqlNU/fastBJM")
library("fastBJM")
```

##  An example
This example illustrates the use of fastBJM to analyse a simulated dataset with the longitudinal and event history data from 5000 people. The multistate submodel is given by the state diagram below. There are five health states (1 - 5) and eight permissible transitions $\Delta=\{(12), (13), (15), (24), (35), (34), (35), (45)\}$. For each person, the first longitudinal measurement is taken at the start of the follow-up and subsequent measurements are taken every 2.5 years until the person exits the study due to either (non-informative) drop-out or death (which is State 5 in the multistate model).

<img width="601" height="467" alt="Screenshot 2026-01-08 at 15 17 29" src="https://github.com/user-attachments/assets/48cfe589-99f5-4913-9102-b1ef2639257a" />

The example code below fits a joint model in which a linear growth curve model with random individual-specific intercepts and slopes is used for the longitudinal submodel. For the multistate submodel, the transition intensities are specified as functions of the effects from two (fixed-effect) covariates, $x_1$ (continuous) and $x_2$ (binary) and a linear association structure with the current value of the underlying longitudinal process. The intensity of a generic permissible transition from State $j$ to State $k$ for person $i$ takes the following structure

$$
\lambda^{(jk)}_i(t)=\lambda^{(jk)}_0(t)\cdot e^{\boldsymbol{X}_i\boldsymbol{\beta}^{(jk)} + a_1(u_i+d_i t)}
$$

where $\boldsymbol{X}_i$ denotes the covariates and $u_i$ and $d_i$ are the intercept and slope. The code below specifies a piecewise constant model on the baseline intensity $\lambda^{(jk)}_0(t)$ with two age-intervals 50-70 and 70-150.

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
##   other settings - leave as default
model_spec$nQs <- 15
model_spec$lmm_correlated_random_effects <- TRUE
model_spec$standardise_y <- T
model_spec$standardise_cnt_x <- T
model_spec$weibull_baseline <- FALSE
model_spec$joint <- !is.null(model_spec$which_longvar)
model_spec$include_fixed_effects <- !is.null(model_spec$X_spec)
if (!model_spec$include_fixed_effects) model_spec$beta_by_transitions <- NULL
model_spec$update_RE_byperson <- TRUE #  update each person independently at the MH step
model_spec$byperson <- TRUE

model_spec$PID <- 'mergeid'  #  person ID that links the event history data to the longitudinal data

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

