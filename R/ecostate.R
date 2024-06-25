#' @title 
#' fit EcoState model 
#'
#' @description 
#' Estimate parameters for an EcoState model
#'
#' @param taxa Character vector of taxa included in model. 
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
#' @param V numeric-matrix with rownames and colnamesmatching \code{taxa}, 
#'        where each element gives the vulnerability parameter for a given
#'        interaction.
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
#' @importFrom TMB config
#' @importFrom checkmate assertDouble
#'
#' @details
#' All \code{taxa} must be included in \code{QB}, \code{PB}, \code{B}, and \code{DC},
#' but additional taxa can be in \code{QB}, \code{PB}, \code{B}, and \code{DC} that
#' are not in \code{taxa}.  So \code{taxa} can be used to redefine the set of modeled
#' species without changing other inputs
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
          EE,
          V,
          fit_B = vector(),
          fit_Q = vector(),
          fit_B0 = vector(),
          fit_EE = vector(),
          fit_eps = vector(),
          control = ecostate_control() ){

  # 
  start_time = Sys.time()
  if( !all(c(fit_B,fit_Q,fit_B0,fit_eps,fit_EE) %in% taxa) ){
    if(isFALSE(control$silent)) warning("Some `fit_B`, `fit_Q`, `fit_B0`, or `fit_eps` not in `taxa`")
  }
  if( any(biomass$Mass==0) ) stop("`biomass$Mass` cannot include zeros, given the assumed lognormal distribution")
  if( any(catch$Mass==0) ) stop("`catch$Mass` cannot include zeros, given the assumed lognormal distribution")

  # Set tmbad.sparse_hessian_compress
  config( tmbad.sparse_hessian_compress = control$tmbad.sparse_hessian_compress, DLL="RTMB" )

  # 
  n_species = length(taxa)
  
  # Configuring inputs
  if(!all(taxa %in% names(PB))) stop("Check names for `PB`")
  if(!all(taxa %in% names(QB))) stop("Check names for `QB`")
  if(!all(taxa %in% names(B))) stop("Check names for `B`")
  if(!all(taxa %in% names(EE))) stop("Check names for `EE`")
  if(!all(taxa %in% rownames(V_ij)) | !all(taxa %in% colnames(V_ij))) stop("Check dimnames for `V`")
  if(!all(taxa %in% rownames(V_ij)) | !all(taxa %in% colnames(V_ij))) stop("Check dimnames for `V`")
  logPB_i = log(PB[taxa])
  logQB_i = log(QB[taxa])
  logB_i = log(B[taxa])
  DC_ij = DC[taxa,taxa]
  EE_i = EE[taxa]
  #noB_i = rep(0,n_species)
  # Indicators 
  which_primary = which( colSums(DC_ij)==0 )
  noB_i = ifelse( is.na(logB_i), 1, 0 )
  
  if(any(is.na(c(logPB_i,logQB_i[-which_primary],DC_ij)))){
    stop("Check `PB` `QB` and `DC` for NAs or `taxa` that are not provided")
  }
  
  if(missing(V)){
    V_ij = array(2, dim=c(n_species,n_species), dimnames=list(taxa,taxa))
    V_ij[,which_primary] = 91  # Default high value from Gaichas et al. 2011
  }else{
    V_ij = V[taxa,taxa]
  }
  Vprime_ij = log(V_ij - 1)
  # V = exp(Vprime) + 1 so 1 < V < Inf
  
  # Rescale DC_ij to sum to 1 by predator
  if( any(abs(colSums(DC_ij)-1) > 0.01) ){
    if(isFALSE(control$silent)) message("Rescaling columns of `DC` to sum to one")
  }
  colsums = colSums(DC_ij)
  DC_ij = DC_ij / outer( rep(1,nrow(DC_ij)), ifelse(colsums==0,1,colsums) )
  
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
  p = list( delta_i = rep(log(1), n_species),
            ln_sdB = log(0.1), 
            ln_sdC = log(0.1),
            logB_i = logB_i,
            EE_i = EE_i,
            logPB_i = logPB_i,
            logQB_i = logQB_i,
            Vprime_ij = Vprime_ij,
            DC_ij = DC_ij,
            logtau_i = rep(NA, n_species),
            epsilon_ti = array( 0, dim=c(nrow(Bobs_ti),n_species) ),
            logF_ti = array( log(0.01), dim=c(nrow(Bobs_ti),n_species) ),
            logq_i = rep( log(1), n_species) )      # , PB_i=PB_i

  # 
  map = list()
  
  # 
  map$logPB_i = factor( rep(NA,n_species) )
  map$logQB_i = factor( rep(NA,n_species) )
  map$DC_ij = factor( array(NA,dim=dim(p$DC_ij)) )
  map$Vprime_ij = factor( array(NA,dim=dim(p$Vprime_ij)) )
  
  # 
  #p$logtau_i = ifelse(taxa %in% fit_eps, log(0.01)+logB_i, NA)
  p$logtau_i = ifelse(taxa %in% fit_eps, log(0.001), NA)
  
  # Catches
  map$logF_ti = factor( ifelse(is.na(Cobs_ti), NA, seq_len(prod(dim(Cobs_ti)))) )
  p$logF_ti[] = ifelse(is.na(Cobs_ti), -Inf, log(0.01))
  
  # Catchability
  map$logq_i = factor( ifelse(taxa %in% fit_Q, seq_len(n_species), NA) )
  
  # Initial biomass-ratio ... turn off if no early biomass observations
  map$delta_i = factor( ifelse(taxa %in% fit_B0, seq_along(p$delta_i), NA) )
  
  # process errors
  map$logtau_i = seq_len(n_species)
  map$epsilon_ti = array( seq_len(prod(dim(p$epsilon_ti))), dim=dim(p$epsilon_ti))
  for(i in 1:n_species){
    if( is.na(p$logtau_i[i]) ){
      map$logtau_i[i] = NA
      map$epsilon_ti[,i] = NA
      p$epsilon_ti[,i] = 0
    }
  }
  map$logtau_i = factor(map$logtau_i)
  map$epsilon_ti = factor(map$epsilon_ti)
  
  # Measurement errors
  p$ln_sdB = log(0.1)
  map$ln_sdB = factor(NA)
  p$ln_sdC = log(0.1)
  map$ln_sdC = factor(NA)

  # Fix biomass for primary producers .... seems to be stiff if trying to fix more than one variable
  map$logB_i = factor( ifelse(taxa %in% fit_B, seq_len(n_species), NA) )
  map$EE_i = factor( ifelse(taxa %in% fit_EE, seq_len(n_species), NA) )
  
  # User-supplied parameters
  if( !is.null(control$tmb_par) ){
    # Check shape but not numeric values, and give informative error
    attr(p,"check.passed") = attr(control$tmb_par,"check.passed")
    if( isTRUE(all.equal(control$tmb_par, p, tolerance=Inf, check.class=FALSE, check.attributes=FALSE)) ){
      message("Using `control$tmb_par`, so be cautious in constructing it")
      p = control$tmb_par
    }else{
      stop("Not using `control$tmb_par` because it has some difference from `p` built internally")
    }
  }
  
  #
  if( !is.null(control$map) ){
    message("Using `control$map`, so be cautious in constructing it")
    map = control$map
  }
  
  # Load data in environment for function "compute_nll"
  data = local({
                  Bobs_ti = Bobs_ti
                  Cobs_ti = Cobs_ti
                  n_steps = control$n_steps 
                  if( control$integration_method == "ABM"){
                    project_vars = abm3pc_sys 
                  }else if( control$integration_method =="RK4"){
                    project_vars = rk4sys 
                  }else{
                    project_vars = function(f, a, b, y0, n, Pars){
                      myode( f, a, b, y0, n, Pars, method=control$integration_method )
                    }
                  }
                  F_type = control$F_type
                  n_species = n_species
                  noB_i = noB_i
                  scale_solver = control$scale_solver
                  inverse_method = control$inverse_method
                  which_primary = which_primary
                  environment()
  })
  environment(compute_nll) <- data

  # Load data in environment for function "dBdt"
  data2 = local({
                  #logPB_i = logPB_i
                  #logQB_i = logQB_i
                  #logB_i = logB_i
                  #DC_ij = DC_ij
                  #logV_ij = logV_ij
                  which_primary = which_primary
                  #EE_i = EE_i
                  n_species = n_species
                  F_type = control$F_type
                  environment()
  })
  environment(dBdt) <- data2
  
  # Make TMB object
  #browser()
  # environment(compute_nll) <- data
  obj <- MakeADFun( func = compute_nll, 
                    parameters = p,
                    map = map,
                    random = c("epsilon_ti"),
                    profile = control$profile,
                    silent = control$silent )
  
  # Optimize
  opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isFALSE(control$quiet) ) message("Running nlminb_loop #", i)
    opt = nlminb( start = opt$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isFALSE(control$quiet) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(opt$par) )
    h = optimHess(opt$par, fn=obj$fn, gr=obj$gr)
    opt$par = opt$par - solve(h, g)
    opt$objective = obj$fn(opt$par)
  }
  rep = obj$report()
  parhat = obj$env$parList()

  # sdreport
  if( isTRUE(control$getsd) ){
    hessian.fixed = optimHess( par = opt$par, 
                      fn = obj$fn, 
                      gr = obj$gr )
    sdrep = sdreport( obj,
                      par.fixed = opt$par,
                      hessian.fixed = hessian.fixed,
                      getJointPrecision = control$getJointPrecision )
  }else{
    hessian.fixed = sdrep = NULL
  }

  #
  environment()
  on.exit( gc() )  # Seems necessary after environment()
  
  # bundle and return output
  internal = list(
    parhat = parhat,
    control = control,
    Bobs_ti = Bobs_ti,
    Cobs_ti = Cobs_ti,
    # Avoid stuff that's in parhat
    #logPB_i = logPB_i, 
    #logQB_i = logQB_i, 
    #logB_i = logB_i, 
    #DC_ij = DC_ij, 
    #Vprime_ij = Vprime_ij,
    hessian.fixed = hessian.fixed,
    taxa = taxa,
    years = years
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
#' @param map list of mapping values, passed to [RTMB::MakeADFun]
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
#' @param integration_method What numerical integration method to use. \code{"ABM"}
#'        uses a native-R versions of Adam-Bashford, code{"RK4"} uses a native-R
#'        version of Runge-Kutta-4, and code{"ode23"} uses a native-R
#'        version of adaptive Runge-Kutta-23, 
#'        where all are adapted from \code{pracma} functions.
#'        \code{"rk4"} and \code{lsoda} use those methods
#'        from \code{deSolve::ode} as implemented by \code{RTMBode::ode}
#' @param F_type whether to integrate catches along with biomass (\code{"integrated"})
#'        or calculate catches from the Baranov catch equation applied to average 
#'        biomass (\code{"averaged"})
#' @param tmbad.sparse_hessian_compress passed to [TMB::config()], and enabling 
#'        an experimental feature to save memory when first computing the inner
#'        Hessian matrix.  Using \code{tmbad.sparse_hessian_compress=1} seems
#'        to have no effect on the MLE (although users should probably confirm this), 
#'        and hugely reduces memory use in both small
#'        and large models. Using \code{tmbad.sparse_hessian_compress=1} seems
#'        to hugely speed up the model-fitting with a large model but results in a small
#'        decrease in speed for model-fitting with a small model. 
#'
#' @export
ecostate_control <-
function( nlminb_loops = 1,
          newton_loops = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          silent = getOption("ecostate.silent", TRUE),
          trace = getOption("ecostate.trace", 1),
          verbose = getOption("ecostate.verbose", FALSE),
          profile = c("logF_ti"),
          tmb_par = NULL,
          map = NULL,
          getJointPrecision = FALSE,
          integration_method = c( "ABM", "RK4", "ode23", "rk4", "lsoda" ),
          n_steps = 10,
          F_type = c("integrated", "averaged"),
          scale_solver = c("joint", "simple"),
          inverse_method = c("Standard", "Penrose_moore"),
          tmbad.sparse_hessian_compress = 1 ){

  #
  integration_method = match.arg(integration_method)
  F_type = match.arg(F_type)
  scale_solver = match.arg(scale_solver)
  inverse_method = match.arg(inverse_method)
  
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
    map = map,
    getJointPrecision = getJointPrecision,
    integration_method = integration_method,
    n_steps = n_steps,
    F_type = F_type,
    scale_solver = scale_solver,
    inverse_method = inverse_method,
    tmbad.sparse_hessian_compress = tmbad.sparse_hessian_compress
  ), class = "ecostate_control" )
}

#' @title Print EcoSim parameters
#'
#' @description Prints parameters defining EcoSim dynamics
#'
#' @param x Output from \code{\link{ecostate}}
#' @param silent whether to print to terminal
#'
#' @return
#' invisibly returns table printed
#'
#' @export
print_ecopars <-
function( x,
          silent = FALSE ){
  
  # Params
  out1 = rbind( 
    "QB" = exp(x$internal$parhat[['logQB_i']]),
    "PB" = exp(x$internal$parhat[['logPB_i']]),
    # Use out_initial so it includes add_equilibrium values
    "B" = x$rep$out_initial$B_i,      
    "EE" = x$rep$out_initial$EE_i
  )
  colnames(out1) = names(x$internal$logPB_i)
  
  # Diet
  out2 = x$internal$parhat[["DC_ij"]]
  
  # Diet
  out3 = exp(x$internal$parhat[["Vprime_ij"]]) + 1
  
  # Print to terminal
  if(isFALSE(silent)){
    cat("EcoSim parameters:\n")
    print(out1)
    cat("\nEcoSim diet matrix:\n")
    print(out2)
    cat("\nEcoSim vulnerability matrix:\n")
    print(out3)
  }
  Return = list( "parameters" = out1, 
                 "diet_matrix" = out2, 
                 "vulnerability_matrix" = out3 )
  return(invisible(Return))
}

#' @title Print fitted ecostate object
#'
#' @description Prints output from fitted ecostate model
#'
#' @param x Output from \code{\link{ecostate}}
#' @param ... Not used
#'
#' @return
#' No return value, called to provide clean terminal output when calling fitted
#' object in terminal.
#'
#' @method print ecostate
#' @export
print.ecostate <-
function( x,
          ... ){
  cat("Dynamics integrated using ", x$internal$control$integration_method, " with ", x$internal$control$n_steps, " time-steps")
  cat("\nRun time: " )
  print(x$run_time)
  cat("Negative log-likelihood: " )
  cat( x$opt$objective )
  cat("\n\n")
  
  # Print pars
  print_ecopars( x )
  
  # Print parameters
  if( !is.null(x$sdrep) ){
    cat("\nEstimates: ")
    print(x$sdrep)
  }
}
