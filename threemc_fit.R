#' Wrapper function for `memprof::with_monitor()` to only return memory_use
#' TMB object for GHA is too large to pull from cluster into local session
memory_use <- function(x) {
  return(x$memory_use)
}

#' Compile and load model from character string
#' (from https://github.com/jeffeaton/inla-sandbox/blob/master/inla-tmb-constraints.R#L12-L24)
#' @param code model code as a character string
#'
#' @return name of the loaded DLL
#'
tmb_compile_and_load <- function(mod, framework = "TMBad", ...) {
  f <- tempfile(fileext = ".cpp", ...)
  if (!dir.exists(dirname(f))) dir.create(dirname(f))
  writeLines(mod, f)
  # for Windows, need to replace "\\" with *Nix-like "/"
  if(grepl("\\", "src\\threemc.cpp", fixed = TRUE)) {
    f <- gsub("\\\\", "/", f)
  }
  TMB::compile(f, framework = framework)
  dyn.load(TMB::dynlib(tools::file_path_sans_ext(f)))
  basename(tools::file_path_sans_ext(f))
}

#' Run example threemc fit for testing
#'
#' @return fit object
#' @export
threemc_fit <- function(shell_dat, areas, mod,  silent = FALSE, ...) {

  # ensure only shapefiles for country in shell_dat are present
  areas <- subset(areas, iso3 %in% substr(shell_dat$area_id, 0, 3))
  # add space column
  areas$space <- seq_len(nrow(areas))

  #### Create model matrices ####

  dat_tmb <- threemc::threemc_prepare_model_data(out = shell_dat, areas = areas)

  #### Create parameters ####

  parameters <- threemc::threemc_initial_pars(dat_tmb = dat_tmb)

  #### Fit TMB model ####

  # compile and load threemc TMB model
  if (length(mod) == 1) {
    dll <- mod
  } else {
    dll <- tmb_compile_and_load(mod, ...)
  }

  # TMB config options
  TMB::config(
    # should reduce memory usage https://tinyurl.com/5cuxmm4t
    # Note: https://github.com/mrc-ide/threemc_example/issues/4
    # best to turn this off, ignored for CppAD & slows down TMBad 
    # tmbad.sparse_hessian_compress = 1,
    # Reduce memory peak of a parallel model by creating tapes in serial
    tape.parallel = 0,
    DLL = dll
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
    DLL        = dll,
    silent     = silent
  )
  
  if (silent == TRUE) {
    obj$env$tracemgc <- FALSE
    obj$env$inner.control$trace <- FALSE
  }

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
