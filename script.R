#### TODO: Add Description Here! ####

#### Libs ####

library(readr)
library(sf)
# devtools::install_github("mrc-ide/threemc")
library(threemc)
library(TMB)

#### Load Data ####

# load shell dataset
shell_dat <- read_csv("data/shell_data.csv.gz")
# load shapefiles
areas <- read_sf("data/areas.geojson")

#### Create model matrices ####

dat_tmb <- threemc_prepare_model_data(out = shell_dat, areas = areas)

#### Create parameters ####

parameters <- threemc_initial_pars(dat_tmb = dat_tmb)

#### Fit TMB model ####

# compile and load threemc TMB model
mod <- "threemc"
compile(paste0(mod, ".cpp"))
dyn.load(dynlib(mod))

# TMB config options 
TMB::config(
  # should reduce memory usage https://tinyurl.com/5cuxmm4t
  tmbad.sparse_hessian_compress = 1, 
  # Reduce memory peak of a parallel model by creating tapes in serial
  tape.parallel = 0, 
  DLL = mod
)

# construct objective function
obj <- MakeADFun(
  data       = dat_tmb, 
  parameters = parameters,  
  random     = c(
    "u_time_mmc", "u_age_mmc", "u_space_mmc",
    "u_agetime_mmc", "u_agespace_mmc", "u_spacetime_mmc",
    "u_age_tmc", "u_space_tmc", 
    "u_agespace_tmc"
  ),
  DLL        = mod
)

# run optimiser (very memory intensive in "optimising tape ..." stage!)
opt <- stats::nlminb(
  start   = obj$par,
  obj     = obj$fn,
  gr      = obj$gr,
  control = list(trace = 1)
)

# sample from fit
fit <- threemc:::circ_sample_tmb(
  obj = obj, opt = opt, nsample = 1000, sdreport = FALSE
)
