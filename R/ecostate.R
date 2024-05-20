#' @title 
#' fit EcoState model 
#'
#' @description 
#' Estimate parameters for an EcoState model
#'
#' @param taxa Character vector of taxa included in model
#' @param years Integer-vector of years included in model                  
#' @param catch long-form data frame with columns \code{Mass}, \code{Year}
#'        and  \code{Taxa}
#' @param biomass long-form data frame with columns \code{Mass}, \code{Year}
#'        and  \code{Taxa}, where \code{Mass} is assumed to have the same units
#'        as \code{catch}
#' @param PB numeric-vector with names matching \code{taxa}, providing the
#'        ratio of production to biomass for each taxon
#' @param QB numeric-vector with names matching \code{taxa}, providing the
#'        ratio of consumption to biomass for each taxon
#' @param B numeric-vector with names matching \code{taxa}, providing the
#'        starting (or fixed) value for equilibrium biomass for each taxon
#' @param DC numeric-matrix with rownames and colnamesmatching \code{taxa}, 
#'        where each column provides the diet proportion for a given predator
#' @param fit_B Character-vector listing \code{taxa} for which equilibrium
#'        biomass is estimated as a fixed effect
#' @param fit_Q Character-vector listing \code{taxa} for which the catchability
#'        coefficient is estimated as a fixed effect
#' @param fit_B0 Character-vector listing \code{taxa} for which the ratio of initial
#'        to equilibrium biomass is estimated as a fixed effect
#' @param fit_eps Character-vector listing \code{taxa} for which the
#'        model should estimate annual process errors in dB/dt
#' @param control Output from [ecostate_control()], used to define user
#'        settings.
#'
#' @details
#' todo
#'
#' @export
ecostate <-
function( taxa,
          years,
          catch,
          biomass,
          PB,
          QB,
          B,
          DC,
          fit_B,
          fit_Q,
          fit_B0,
          fit_eps,
          control = ecostate_control() ){

  # 
  n_species = length(taxa)
  
  # Configuring inputs
  logPB_i = log(PB[taxa])
  logQB_i = log(QB[taxa])
  logB_i = log(B[taxa])
  DC_ij = DC[taxa,taxa]
  logV_ij = matrix( log(2), nrow=n_species, ncol=n_species )
  
  # Convert long-form `catch` to wide-form Cobs_ti
  Cobs_ti = tapply( catch[,'Mass'], FUN=mean, INDEX = list(
                    factor(catch[,'Year'],levels=years),
                    factor(catch[,'Taxon'],levels=taxa) )
                  )
  if(any(!is.na(Cobs_ti[1,]))) message("Fixing catch=NA in first year as required")
  Cobs_ti[1,] = NA
  
  # Convert long-form `biomass` to wide-form Bobs_ti
  Bobs_ti = tapply( biomass[,'Mass'], FUN=mean, INDEX = list(
                    factor(biomass[,'Year'],levels=years),
                    factor(biomass[,'Taxon'],levels=taxa) )
                  )
  
  # parameter list
  p = list( logB0ratio_i = rep(log(1), n_species),
            ln_sdB = log(0.1), 
            ln_sdC = log(0.1),
            logB_i = log(B0_i),
            logtau_i = rep(NA, n_species),
            deltaB_ti = rarray(Bobs_ti, sd=0),
            logF_ti = array( log(0.01), dim=c(nrow(Bobs_ti),n_species) ),
            logq_i = rep( log(1), n_species) )      # , PB_i=PB_i

  # 
  map = list()
  
  # 
  p$logtau_i = ifelse(taxa %in% fit_delta, log(0.01*B0_i), NA)
  
  # Catches
  map$logF_ti = factor( ifelse(is.na(Cobs_ti), NA, seq_len(prod(dim(Cobs_ti)))) )
  p$logF_ti[] = ifelse(is.na(Cobs_ti), -Inf, log(0.01))
  
  # Catchability
  map$logq_i = factor( ifelse(taxa %in% fit_Q, seq_len(n_species), NA) )
  
  # Initial biomass-ratio ... turn off if no early biomass observations
  map$logB0ratio_i = factor( ifelse(taxa %in% fit_B0, seq_along(p$logB0ratio_i), NA) )
  
  # process errors
  map$logtau_i = seq_len(n_species)
  map$deltaB_ti = array( seq_len(prod(dim(p$deltaB_ti))), dim=dim(p$deltaB_ti))
  for(i in 1:n_species){
    if( is.na(p$logtau_i[i]) ){
      map$logtau_i[i] = NA
      map$deltaB_ti[,i] = NA
      p$deltaB_ti[,i] = 0
    }
  }
  map$logtau_i = factor(map$logtau_i)
  map$deltaB_ti = factor(map$deltaB_ti)
  
  # Measurement errors
  p$ln_sdB = log(0.1)
  map$ln_sdB = factor(NA)
  p$ln_sdC = log(0.1)
  map$ln_sdC = factor(NA)
  
  # Fix biomass for primary producers .... seems to be stiff if trying to fix more than one variable
  map$logB_i = factor( ifelse(taxa %in% fit_B, seq_len(n_species), NA) )
  
  # Make TMB object
  #browser()
  obj <- MakeADFun( func = compute_nll, 
                    parameters = p,
                    map = map,
                    random = c("deltaB_ti"),
                    profile = c("logF_ti"),
                    silent = TRUE )
  #obj$fn(obj$par); obj$gr(obj$par)
  
  # 
  start_time = Sys.time()
  opt = nlminb( start = obj$par, 
                obj = obj$fn, 
                gr = obj$gr,
                control = list(eval.max=control$eval.max, iter.max=control$iter.max, trace=control$trace) )
  run_time = Sys.time() - start_time
  sdrep = sdreport(obj)              
  rep = obj$report()
  parhat = obj$env$parList()

  # bundle and return output
  internal = list(
    parhat = parhat,
    control = control
  )
  out = list(
    obj = obj,
    opt = opt,
    rep = rep,
    sdrep = sdrep,
    tmb_inputs = list(p=p, map=map),
    call = match.call(),
    run_time = Sys.time() - start_time,
    internal = internal
  )

  class(out) = "ecostate"
  return(out)
} 

#' @title Detailed control for ecostate structure
#'
#' @description Define a list of control parameters.  
#'
#' @param nlminb_loops Integer number of times to call [stats::nlminb()].
#' @param newton_loops Integer number of Newton steps to do after running
#'   [stats::nlminb()].
#' @param getsd Boolean indicating whether to call [TMB::sdreport()]
#' @param tmb_par list of parameters for starting values, with shape identical
#'   to `tinyVAST(...)$internal$parlist`
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed. Passed to `control` in [stats::nlminb()].
#' @param iter.max Maximum number of iterations allowed. Passed to `control` in
#'   [stats::nlminb()].
#' @param verbose Output additional messages about model steps during fitting?
#' @param silent Disable terminal output for inner optimizer?
#' @param trace Parameter values are printed every `trace` iteration
#'   for the outer optimizer. Passed to
#'   `control` in [stats::nlminb()].
#' @param getJointPrecision whether to get the joint precision matrix.  Passed
#'        to \code{\link[TMB]{sdreport}}.
#'
#' @export
ecostate_control <-
function( nlminb_loops = 1,
          newton_loops = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          silent = getOption("ecostate.silent", TRUE),
          trace = getOption("ecostate.trace", 0),
          verbose = getOption("ecostate.verbose", FALSE),
          profile = c(),
          tmb_par = NULL,
          getJointPrecision = FALSE ){

  # Return
  structure( list(
    nlminb_loops = nlminb_loops,
    newton_loops = newton_loops,
    eval.max = eval.max,
    iter.max = iter.max,
    getsd = getsd,
    silent = silent,
    trace = trace,
    verbose = verbose,
    profile = profile,
    tmb_par = tmb_par,
    getJointPrecision = getJointPrecision
  ), class = "ecostate_control" )
}

