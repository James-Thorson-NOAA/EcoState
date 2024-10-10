#' @title 
#' fit EcoState model 
#'
#' @description 
#' Estimate parameters for an EcoState model
#'
#' @param taxa Character vector of taxa included in model. 
#' @param years Integer-vector of years included in model                  
#' @param catch long-form data frame with columns \code{Mass}, \code{Year}
#'        and  \code{Taxon}
#' @param biomass long-form data frame with columns \code{Mass}, \code{Year}
#'        and  \code{Taxon}, where \code{Mass} is assumed to have the same units
#'        as \code{catch}
#' @param PB numeric-vector with names matching \code{taxa}, providing the                        
#'        ratio of production to biomass for each taxon
#' @param QB numeric-vector with names matching \code{taxa}, providing the          
#'        ratio of consumption to biomass for each taxon
#' @param B numeric-vector with names matching \code{taxa}, providing the
#'        starting (or fixed) value for equilibrium biomass for each taxon
#' @param U numeric-vector with names matching \code{taxa}, providing the 
#'        proportion of consumption that is unassimilated and therefore
#'        exported to detritus
#' @param type character-vector with names matching \code{taxa} and
#'        elements \code{c("auto","hetero","detritus")},
#'        indicating whether each taxon is a primary producer, consumer/predator, or
#'        detritus, respectively.   
#' @param DC numeric-matrix with rownames and colnames matching \code{taxa}, 
#'        where each column provides the diet proportion for a given predator
#' @param X numeric-matrix with rownames and colnames matching \code{taxa}, 
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
#' @param fit_PB Character-vector listing \code{taxa} for which equilibrium
#'        production per biomass is estimated.  Note that it is likely
#'        a good idea to include a prior for any species for which this is estimated.
#' @param fit_QB Character-vector listing \code{taxa} for which equilibrium
#'        consumption per biomass is estimated.  Note that it is likely
#'        a good idea to include a prior for any species for which this is estimated.
#' @param log_prior A user-provided function that takes as input the list of
#'        parameters \code{out$obj$env$parList()} where \code{out} is the output from
#'        \code{ecostate()}, and returns a numeric vector
#'        where the sum is the log-prior probability.  For example
#'        \code{log_prior = function(p) dnorm( p$logq_i[1], mean=0, sd=0.1, log=TRUE)}
#'        specifies a lognormal prior probability for the catchability coefficient
#'        for the first \code{taxa} with logmean of zero and logsd of 0.1
#' @param control Output from [ecostate_control()], used to define user
#'        settings.
#'
#' @importFrom TMB config
#' @importFrom checkmate assertDouble assertFactor assertCharacter assertList
#' @importFrom Matrix Matrix Diagonal sparseMatrix
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
          catch = data.frame("Year"=numeric(0),"Mass"=numeric(0),"Taxon"=numeric(0)),
          biomass = data.frame("Year"=numeric(0),"Mass"=numeric(0),"Taxon"=numeric(0)),
          agecomp = list(),
          weight = list(),
          PB,
          QB,
          B,
          DC,
          EE,
          X,
          type,
          U,
          fit_B = vector(),
          fit_Q = vector(),
          fit_B0 = vector(),
          fit_EE = vector(),
          fit_PB = vector(),
          fit_QB = vector(),
          fit_eps = vector(),
          fit_nu = vector(),
          fit_phi = vector(),
          log_prior = function(p) 0,
          settings = stanza_settings(taxa=taxa),
          control = ecostate_control() ){

  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

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
  
  if(missing(U)){
    U = rep(0.2, n_species)
    names(U) = taxa
  }  
  if(missing(type)){
    type = ifelse(colSums(DC_ij)==0, "auto", "hetero")
    names(type) = taxa
  }
  
  # Configuring inputs
  if(!all(taxa %in% names(PB))) stop("Check names for `PB`")
  if(!all(taxa %in% names(QB))) stop("Check names for `QB`")
  if(!all(taxa %in% names(B))) stop("Check names for `B`")
  if(!all(taxa %in% names(EE))) stop("Check names for `EE`")
  if(!all(taxa %in% names(type))) stop("Check names for `type`")
  if(!all(taxa %in% names(U))) stop("Check names for `U`")
  logPB_i = array( log(PB[taxa]), dimnames=list(taxa) )
  logQB_i = log(QB[taxa])
  logB_i = log(B[taxa])
  DC_ij = DC[taxa,taxa,drop=FALSE]
  EE_i = EE[taxa]
  type_i = type[taxa]
  U_i = U[taxa]
  
  # Deal with V
  if(missing(X)){
    X_ij = array(2, dim=c(n_species,n_species), dimnames=list(taxa,taxa))
    X_ij[,which_primary,drop=FALSE] = 91  # Default high value from Gaichas et al. 2011
  }else{
    if(!all(taxa %in% rownames(X)) | !all(taxa %in% colnames(X))) stop("Check dimnames for `X`")
    if(!all(taxa %in% rownames(X)) | !all(taxa %in% colnames(X))) stop("Check dimnames for `X`")
    X_ij = X[taxa,taxa,drop=FALSE]
  }
  Xprime_ij = log(X_ij - 1)
  # V = exp(Xprime) + 1 so 1 < V < Inf
  
  #
  assertCharacter( type_i, len=n_species, any.missing=FALSE )
  if(isFALSE(all(type_i %in% c("auto","hetero","detritus")))) stop("Confirm ", type, " only contains auto, hetero, or detritus")
  if(sum(type_i=="detritus") >=2) stop("Currently can only specify one detritus variable")
  assertDouble( U_i, len=n_species, any.missing=FALSE, upper=1 )      # GE = 1-U-A and A>=0 so GE <= 1-U so GE+U <= 1

  #noB_i = rep(0,n_species)
  # Indicators 
  which_primary = which( type_i=="auto" )
  which_detritus = which( type_i=="detritus" )
  which_multigroup = match( settings$multigroup_taxa, settings$taxa )
  noB_i = ifelse( is.na(logB_i), 1, 0 )
  noB_i[which_multigroup] = 0

  if(any(is.na(c(logPB_i,logQB_i[-c(which_primary,which_detritus,which_multigroup)],DC_ij)))){
    stop("Check `PB` `QB` and `DC` for NAs or `taxa` that are not provided")
  }
  
  # Rescale DC_ij to sum to 1 by predator
  if( any(abs(colSums(DC_ij)-1) > 0.01) ){
    if(isFALSE(control$silent)) message("Rescaling columns of `DC` to sum to one")
  }
  colsums = colSums(DC_ij)
  DC_ij = DC_ij / outer( rep(1,nrow(DC_ij)), ifelse(colsums==0,1,colsums) )
  
  # Convert to sparse diet matrix
  #DC_ij = Matrix::Matrix(DC_ij)
  
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
  
  #
  assertList( agecomp )
  Nobs_ta_g2 = agecomp[match(names(agecomp),settings$unique_stanza_groups)]  # match works for empty list
  # ADD MORE CHECKS

  #
  assertList( weight )
  Wobs_ta_g2 = weight[match(names(weight),settings$unique_stanza_groups)]  # match works for empty list

  #
  stanza_data = make_stanza_data( settings )

  # number of selex params
  n_selex = length(Nobs_ta_g2) # * switch( settings$comp_weight, "multinom" = 2, "dir" = 3, "dirmult" = 3 )
  n_weight = length(Wobs_ta_g2)

  # parameter list
  #browser()
  p = list( delta_i = rep(log(1), n_species),
            ln_sdB = log(0.1), 
            ln_sdC = log(0.1),
            logB_i = logB_i,
            EE_i = EE_i,
            logPB_i = logPB_i,
            logQB_i = logQB_i,
            U_i = U_i,
            Xprime_ij = Xprime_ij,
            DC_ij = DC_ij,
            logtau_i = rep(NA, n_species),
            logsigma_i = rep(NA, n_species),
            logpsi_g2 = rep(NA, settings$n_g2),
            epsilon_ti = array( 0, dim=c(0,n_species) ),
            alpha_ti = array( 0, dim=c(0,n_species) ),
            nu_ti = array( 0, dim=c(0,n_species) ),
            phi_tg2 = array( 0, dim=c(0,settings$n_g2) ),
            logF_ti = array( log(0.01), dim=c(nrow(Bobs_ti),n_species) ),
            logq_i = rep( log(1), n_species),
            s50_z = rep(1, n_selex),
            srate_z = rep(1, n_selex),
            compweight_z = rep(1, n_selex*ifelse(settings$comp_weight=="multinom",0,1)),
            #selex_z = rep(1, n_selex),  # CHANGE WITH NUMBER OF PARAMETERS
            log_winf_z = rep(0, n_weight),
            ln_sdW_z = rep(0, n_weight),
            SpawnX_g2 = stanza_data$stanzainfo_g2z[,'SpawnX'],
            log_K_g2 = log(stanza_data$stanzainfo_g2z[,'K']),
            logit_d_g2 = qlogis(stanza_data$stanzainfo_g2z[,'d']),
            Wmat_g2 = stanza_data$stanzainfo_g2z[,'Wmat']
  )      # , PB_i=PB_i

  # 
  map = list()
  
  # map these off by default
  map$U_i = factor( rep(NA,n_species) )
  map$DC_ij = factor( array(NA,dim=dim(p$DC_ij)) )
  map$Xprime_ij = factor( array(NA,dim=dim(p$Xprime_ij)) )
  map$SpawnX_g2 = factor( rep(NA,length(p$SpawnX_g2)) )
  map$Wmat_g2 = factor( rep(NA,length(p$Wmat_g2)) )

  # map these off based on options
  map$logQB_i = factor( ifelse(taxa %in% fit_QB, seq_len(n_species), NA) )
  map$logPB_i = factor( ifelse(taxa %in% fit_PB, seq_len(n_species), NA) )
  map$log_K_g2 = factor(ifelse(settings$unique_stanza_groups %in% settings$fit_K, seq_len(settings$n_g2), NA)) #  factor( rep(NA,length(p$K_g2)) )
  map$logit_d_g2 = factor(ifelse(settings$unique_stanza_groups %in% settings$fit_d, seq_len(settings$n_g2), NA)) #  factor( rep(NA,length(p$K_g2)) )

  #
  p$logtau_i = ifelse(taxa %in% fit_eps, log(control$start_tau), NA)
  map$logtau_i = factor(ifelse(taxa %in% fit_eps, seq_len(n_species), NA))
  
  #
  p$logsigma_i = ifelse(taxa %in% fit_nu, log(control$start_tau), NA)
  map$logsigma_i = factor(ifelse(taxa %in% fit_nu, seq_len(n_species), NA))

  #
  p$logpsi_g2 = ifelse(settings$unique_stanza_groups %in% fit_phi, log(control$start_tau), NA)
  map$logpsi_g2 = factor(ifelse(settings$unique_stanza_groups %in% fit_phi, seq_len(settings$n_g2), NA))

  # Catches
  map$logF_ti = factor( ifelse(is.na(Cobs_ti), NA, seq_len(prod(dim(Cobs_ti)))) )
  p$logF_ti[] = ifelse(is.na(Cobs_ti), -Inf, log(0.01))
  
  # Catchability
  map$logq_i = factor( ifelse(taxa %in% fit_Q, seq_len(n_species), NA) )
  
  # Initial biomass-ratio ... turn off if no early biomass observations
  map$delta_i = factor( ifelse(taxa %in% fit_B0, seq_along(p$delta_i), NA) )
  
  # process errors
  if( control$process_error == "epsilon" ){
    p$epsilon_ti = array( 0, dim=c(nrow(Bobs_ti),n_species) )
    map$epsilon_ti = array( seq_len(prod(dim(p$epsilon_ti))), dim=dim(p$epsilon_ti))
    for(i in seq_len(n_species)){
      if( is.na(p$logtau_i[i]) ){
        p$epsilon_ti[,i] = 0
        map$epsilon_ti[,i] = NA
      }
    }
    map$epsilon_ti = factor(map$epsilon_ti)
  }else{
    p$alpha_ti = array( 0, dim=c(nrow(Bobs_ti),n_species) )
    map$alpha_ti = array( seq_len(prod(dim(p$alpha_ti))), dim=dim(p$alpha_ti))
    for(i in seq_len(n_species)){
      if( is.na(p$logtau_i[i]) ){
        p$alpha_ti[,i] = 0
        map$alpha_ti[,i] = NA
      }
    }
    map$alpha_ti = factor(map$alpha_ti)
  }
  # Variation in consumption
  p$nu_ti = array( 0, dim=c(nrow(Bobs_ti),n_species) )
  map$nu_ti = array( seq_len(prod(dim(p$nu_ti))), dim=dim(p$nu_ti))
  for(i in seq_len(n_species)){
    if( is.na(p$logsigma_i[i]) ){
      p$nu_ti[,i] = 0
      map$nu_ti[,i] = NA
    }
  }
  map$nu_ti = factor(map$nu_ti)
  # Variation in recruitment
  p$phi_tg2 = array( 0, dim=c(nrow(Bobs_ti),settings$n_g2) )
  map$phi_tg2 = array( seq_len(prod(dim(p$phi_tg2))), dim=dim(p$phi_tg2))
  for(g2 in seq_len(settings$n_g2)){
    if( is.na(p$logpsi_g2[g2]) ){
      p$phi_tg2[,g2] = 0
      map$phi_tg2[,g2] = NA
    }
  }
  map$phi_tg2 = factor(map$phi_tg2)

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
                  Nobs_ta_g2 = Nobs_ta_g2
                  Wobs_ta_g2 = Wobs_ta_g2
                  years = years
                  #DC_ij = DC_ij
                  control = control
                  #n_steps = control$n_steps
                  if( control$integration_method == "ABM"){
                    project_vars = abm3pc_sys 
                  }else if( control$integration_method =="RK4"){
                    project_vars = rk4sys 
                  }else{
                    project_vars = function(f, a, b, y0, n, Pars){
                      myode( f, a, b, y0, n, Pars, method=control$integration_method )
                    }
                  }
                  #F_type = control$F_type
                  n_species = n_species
                  noB_i = noB_i
                  #scale_solver = control$scale_solver
                  #inverse_method = control$inverse_method
                  type_i = type_i
                  #process_error = control$process_error
                  #sdreport_detail = control$sdreport_detail
                  settings = settings
                  stanza_data = stanza_data
                  taxa = taxa
                  fit_eps = fit_eps
                  fit_nu = fit_nu
                  log_prior = log_prior
                  environment()
  })
  environment(compute_nll) <- data

  # Load data in environment for function "dBdt"
  data2 = local({
                  type_i = type_i
                  n_species = n_species
                  F_type = control$F_type
                  environment()
  })
  environment(dBdt) <- data2
  environment(project_stanzas) <- data2    # project_stanzas(.) calls dBdt(.)

  # 
  #data3 = local({
  #                DC_ij = DC_ij
  #                environment()
  #})
  #environment(add_equilibrium) <- data3

  # Load data in environment for function "dBdt"
  data4 = local({
                  "c" <- ADoverload("c")
                  "[<-" <- ADoverload("[<-")
                  environment()
  })
  environment(log_prior) <- data4

  # Make TMB object
  #browser()
  # compute_nll(p)
  # environment(compute_nll) <- data
  obj <- MakeADFun( func = compute_nll, 
                    parameters = p,
                    map = map,
                    random = control$random,
                    profile = control$profile,
                    silent = control$silent )
  #traceback(max=20)
  #obj$fn(obj$par)

  # Optimize
  opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isFALSE(control$quiet) ) message("Running nlminb_loop #", i)
    opt = nlminb( start = opt$par,
                  objective = obj$fn,
                  gradient = list(obj$gr,NULL)[[ifelse(control$use_gradient,1,2)]],
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }

  #
  get_hessian = function(obj, par){
    if( (length(obj$env$random)==0) & (sum(obj$env$profile)==0) ){
      H = obj$he(x=par)
    }else{
      H = optimHess(par, fn=obj$fn, gr=obj$gr)
    }
    return(H)
  }

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isFALSE(control$quiet) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(opt$par) )
    #h = optimHess(opt$par, fn=obj$fn, gr=obj$gr)
    h = get_hessian(obj=obj, par=opt$par)
    opt$par = opt$par - solve(h, g)
    opt$objective = obj$fn(opt$par)
  }
  rep = obj$report()
  parhat = obj$env$parList()

  # Sanity checks
  if( any(rep$B_ti==0) ){
    warning("Some `B_ti=0` which typically occurs in multistanza models when W<Wmat for one or more years")
  }

  # sdreport
  derived = list()
  if( isTRUE(control$getsd) ){
    #hessian.fixed = optimHess( par = opt$par, 
    #                  fn = obj$fn, 
    #                  gr = obj$gr )
    hessian.fixed = get_hessian(obj=obj, par=opt$par)
    sdrep = sdreport( obj,
                      par.fixed = opt$par,
                      hessian.fixed = hessian.fixed,
                      getJointPrecision = control$getJointPrecision )
    if( length(control$derived_quantities) > 0 ){
      # Relist
      tmp_list = as.list(sdrep, report=TRUE, what="Estimate")[["do.call(\"c\", derived_values)"]]
      derived$Est = relist( flesh=tmp_list, skeleton=obj$report()$derived_values )
      tmp_list = as.list(sdrep, report=TRUE, what="Std. Error")[["do.call(\"c\", derived_values)"]]
      derived$SE = relist( flesh=tmp_list, skeleton=obj$report()$derived_values )
    }
  }else{
    hessian.fixed = sdrep = NULL
  }

  #
  environment()
  on.exit( gc() )  # Seems necessary after environment()
  
  # bundle and return output (all necessary inputs for compute_nll)
  internal = list(
    call = match.call(),
    parhat = parhat,
    control = control,
    settings = settings,
    log_prior = log_prior,
    Bobs_ti = Bobs_ti,
    Cobs_ti = Cobs_ti,
    Nobs_ta_g2 = Nobs_ta_g2,
    Wobs_ta_g2 = Wobs_ta_g2,
    n_species = n_species,
    noB_i = noB_i,
    biomass = biomass,
    catch = catch,
    # Avoid stuff that's in parhat
    #logPB_i = logPB_i, 
    #logQB_i = logQB_i, 
    #logB_i = logB_i, 
    #DC_ij = DC_ij, 
    #Xprime_ij = Xprime_ij,
    hessian.fixed = hessian.fixed,
    taxa = taxa,
    years = years,
    type_i = type_i
  )
  out = list(
    obj = obj,
    opt = opt,
    rep = rep,
    sdrep = sdrep,
    derived = derived,
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
#' @param process_error Whether to include process error as a continuous rate
#'        (i.e., an "innovation" parameterization, \code{process_error="epsilon"}) 
#'        or as a discrete difference between expected
#'        and predicted biomass (i.e., a "state-space" parameterization),  
#'        \code{process_error="alpha"}The
#'        former is more interpretable, whereas the latter is much more computationally
#'        efficient.  
#' @param F_type whether to integrate catches along with biomass (\code{"integrated"})
#'        or calculate catches from the Baranov catch equation applied to average 
#'        biomass (\code{"averaged"})
#' @param derived_quantities character-vector listing objects to ADREPORT
#' @param tmbad.sparse_hessian_compress passed to [TMB::config()], and enabling 
#'        an experimental feature to save memory when first computing the inner
#'        Hessian matrix.  Using \code{tmbad.sparse_hessian_compress=1} seems
#'        to have no effect on the MLE (although users should probably confirm this), 
#'        and hugely reduces memory use in both small
#'        and large models. Using \code{tmbad.sparse_hessian_compress=1} seems
#'        to hugely speed up the model-fitting with a large model but results in a small
#'        decrease in speed for model-fitting with a small model. 
#' @param start_tau Starting value for the standard deviation of process errors
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
          profile = c("logF_ti","log_winf_z","s50_z","srate_z"),
          random = c("epsilon_ti","alpha_ti","nu_ti","phi_tg2"),
          tmb_par = NULL,
          map = NULL,
          getJointPrecision = FALSE,
          integration_method = c( "ABM", "RK4", "ode23", "rk4", "lsoda" ),
          process_error = c("epsilon", "alpha"),
          n_steps = 10,
          F_type = c("integrated", "averaged"),
          derived_quantities = c("h_g2","B_ti","B0_i"),
          scale_solver = c("joint", "simple"),
          inverse_method = c("Standard", "Penrose_moore"),
          tmbad.sparse_hessian_compress = 1,
          use_gradient = TRUE,
          start_tau = 0.001 ){

  #
  integration_method = match.arg(integration_method)
  F_type = match.arg(F_type)
  scale_solver = match.arg(scale_solver)
  inverse_method = match.arg(inverse_method)
  process_error = match.arg(process_error)
  
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
    random = random,
    tmb_par = tmb_par,
    map = map,
    getJointPrecision = getJointPrecision,
    integration_method = integration_method,
    n_steps = n_steps,
    F_type = F_type,
    derived_quantities = derived_quantities,
    scale_solver = scale_solver,
    inverse_method = inverse_method,
    process_error = process_error,
    tmbad.sparse_hessian_compress = tmbad.sparse_hessian_compress,
    use_gradient = use_gradient,
    start_tau = start_tau
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
  out1 = data.frame( 
    "type" = x$internal$type_i,
    # Use out_initial so it includes add_equilibrium values
    "QB" = x$rep$out_initial$QB_i,
    "PB" = x$rep$out_initial$PB_i,
    "B" = x$rep$out_initial$B_i,      
    "EE" = x$rep$out_initial$EE_i,
    "U" = x$internal$parhat[["U_i"]]
  )
  #colnames(out1) = x$internal$taxa
  
  # Diet
  out2 = x$internal$parhat[["DC_ij"]]
  
  # Vulnerability
  out3 = exp(x$internal$parhat[["Xprime_ij"]]) + 1
  
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

#' @title Marginal log-likelihood
#'
#' @description Extract the (marginal) log-likelihood of a ecostate model
#'
#' @param object Output from \code{\link{ecostate}}
#' @param ... Not used
#'
#' @return object of class \code{logLik} with attributes
#'   \item{val}{log-likelihood}
#'   \item{df}{number of parameters}
#' @importFrom stats logLik
#'
#' @return
#' Returns an object of class logLik. This has attributes
#' "df" (degrees of freedom) giving the number of (estimated) fixed effects
#' in the model, abd "val" (value) giving the marginal log-likelihood.
#' This class then allows \code{AIC} to work as expected.
#'
#' @export
logLik.ecostate <- function(object, ...) {
  val = -1 * object$opt$objective
  df = length( object$opt$par )
  out = structure( val,
             df = df,
             class = "logLik")
  return(out)
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
