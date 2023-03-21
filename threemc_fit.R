#' Run example threemc fit for testing
#'
#' @return fit object
#' @export
threemc_fit <- function() {
  # load shell dataset
  shell_dat <- readr::read_csv("data/shell_data.csv.gz")
  # load shapefiles
  areas <- sf::read_sf("data/areas.geojson")

  #### Create model matrices ####

  dat_tmb <- threemc::threemc_prepare_model_data(out = shell_dat, areas = areas)

  #### Create parameters ####

  parameters <- threemc::threemc_initial_pars(dat_tmb = dat_tmb)

  #### Fit TMB model ####

  # compile and load threemc TMB model
  mod <- "threemc"
  TMB::compile(paste0(mod, ".cpp"))
  dyn.load(TMB::dynlib(mod))

  # TMB config options
  TMB::config(
    # should reduce memory usage https://tinyurl.com/5cuxmm4t
    tmbad.sparse_hessian_compress = 1,
    # Reduce memory peak of a parallel model by creating tapes in serial
    tape.parallel = 0,
    DLL = mod
  )

  # construct objective function
  obj <- TMB::MakeADFun(
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
  threemc:::circ_sample_tmb(
    obj = obj, opt = opt, nsample = 1000, sdreport = FALSE
  )
}
