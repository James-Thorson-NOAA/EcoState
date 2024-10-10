#' @title 
#' Compute negative log-likelihood for EcoState model 
#'
#' @description                                               
#' Compute negative log-likelihood for EcoState model
#'
#' @param p list of parameters
#'
#' @details
#' todo
#'
#' @export
compute_nll <-
function( p ) { 
  
  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  # Compute stanza stuff
  p = add_stanza_params( p,
                   stanza_data = stanza_data,
                   settings = settings )

  # Compute equilibrium values
  p = add_equilibrium( p,
                       scale_solver = control$scale_solver,
                       noB_i = noB_i,
                       type_i = type_i )
  
  # Extract epsilon_ti (local copy be modified later)
  epsilon_ti = p$epsilon_ti
  if(control$process_error=="alpha"){
    epsilon_ti = matrix( 0, ncol=n_species, nrow=nrow(Bobs_ti) )
  }
  
  #
  F_type = control$F_type

  # unfished M0 and M2, and B_i solved for EE_i
  #browser()
  p_t = p
    p_t$epsilon_i = rep(0,n_species)
    p_t$logF_i = rep(-Inf,n_species)
    p_t$nu_i = rep(0,n_species)
    p_t$phi_g2 = rep(0,settings$n_g2)
  out_initial = dBdt( Time = 1,
              State = c( exp(p$logB_i), rep(0,n_species)),
              Pars = p_t,
              what = "stuff" )
  
  # Objects to save
  TL_ti = dBdt0_ti = M_ti = m_ti = G_ti = g_ti = M2_ti = m2_ti = Bmean_ti = Chat_ti = B_ti = Bhat_ti = matrix( NA, ncol=n_species, nrow=nrow(Bobs_ti) )
  loglik1_ti = loglik2_ti = loglik3_ti = loglik4_ti = matrix( 0, ncol=n_species, nrow=nrow(Bobs_ti) )  # Missing = 0
  loglik5_tg2 = loglik6_tg2 = loglik7_tg2 = matrix( 0, nrow=nrow(Bobs_ti), ncol=length(settings$unique_stanza_groups) )
  Q_tij = array( NA, dim=c(nrow(Bobs_ti),n_species,n_species) )
  Nexp_ta_g2 = Nobs_ta_g2
  Wexp_ta_g2 = Wobs_ta_g2
  TotalSB_tg2 = TotalEggs_tg2 = matrix( 0, nrow=nrow(Bobs_ti), ncol=stanza_data$n_g2 )
  TotalZ_ts2 = matrix( 0, nrow=nrow(Bobs_ti), ncol=stanza_data$n_s2 )

  # Define initial condition
  B_ti[1,] = out_initial$B_i * exp(p$delta_i)
  TotalEggs_tg2[1,] = p$baseEggs_g2 * settings$STEPS_PER_YEAR
  TotalSB_tg2[1,] = p$baseSB_g2 * settings$STEPS_PER_YEAR
  TotalZ_ts2[1,] = 0
  jnll = 0
  if( control$process_error == "alpha" ){
    epsilon_ti[1,] = p$alpha_ti[1,] 
  }
  Y_zz_g2 = p$Y_zz_g2

  # Extract initial conditions for later reporting
  B0_i = out_initial$B_i
  TL0_i = compute_tracer( Q_ij = out_initial$Q_ij,
                              inverse_method = control$inverse_method,
                              type_i = type_i,
                              tracer_i = rep(1,n_species) )
  G0_ti = out_initial$G_i
  g0_ti = out_initial$g_i
  M0_ti = out_initial$M_i
  m0_ti = out_initial$m_i
  M20_ti = out_initial$M2_i
  m20_ti = out_initial$m2_i

  # 
  Y_tzz_g2 = NULL
  for( g2 in seq_along(Y_zz_g2) ){
    Y_tzz_g2[[g2]] = array( 0.0, dim=c(nrow(Bobs_ti),dim(p$Y_zz_g2[[g2]])),
                 dimnames=list(NULL,rownames(p$Y_zz_g2[[g2]]),colnames(p$Y_zz_g2[[g2]])) )
    Y_tzz_g2[[g2]][1,,] = Y_zz_g2[[g2]]
  }
  #Y_tzz = array( 0.0, dim=c(nrow(Bobs_ti),dim(p$Y_zz)),
  #               dimnames=list(NULL,rownames(p$Y_zz),colnames(p$Y_zz)) ) # dimnames-list dimension matching
  #Y_tzz[1,,] = Y_zz

  # Loop through years
  for( t in 2:nrow(B_ti) ){
  #for( t in 2:66 ){
    # Assemble inputs
    p_t = p
    p_t$Y_zz_g2 = Y_zz_g2
    p_t$logF_i = p$logF_ti[t,]

    # State-space or continuous innovations
    if( control$process_error == "epsilon" ){
      p_t$epsilon_i = epsilon_ti[t,]
    }else{
      p_t$epsilon_i = rep(0,n_species)
    }
    p_t$nu_i = p$nu_ti[t,]
    p_t$phi_g2 = p$phi_tg2[t,]

    # RTMBode::ode requires y0 have names
    y0 = c(B_ti[t-1,], rep(0,n_species))
    names(y0) = paste0("var_",seq_along(y0))

    # Project dynamics
    #browser()
    proj = project_vars(
          f = dBdt,
          a = 0, 
          b = 1,
          n = control$n_steps,
          Pars = p_t,
          y0 = y0 )
    Bnew_i = proj$y[nrow(proj$y),seq_len(n_species)]

    if( settings$n_g2 > 0 ){
      # Project stanzas
      proj_stanzas = project_stanzas(
                   p = p_t,
                   stanza_data = stanza_data,
                   y = proj$y,
                   STEPS_PER_YEAR = settings$STEPS_PER_YEAR,
                   record_steps = FALSE,
                   correct_errors = settings$correct_errors ) # Have used TRUE
      #Y_zz = proj_stanzas$Y_zz
      #Y_tzz[t,,] = Y_zz
      Y_zz_g2 = proj_stanzas$Y_zz_g2
      #browser()
      TotalEggs_tg2[t,] = proj_stanzas$TotalEggs_g2
      TotalSB_tg2[t,] = proj_stanzas$TotalSB_g2
      TotalZ_ts2[t,] = proj_stanzas$TotalZ_s2
      for( g2 in seq_along(Y_zz_g2) ){
        Y_tzz_g2[[g2]][t,,] = Y_zz_g2[[g2]]
      }

      #
      Bnew_s2 = get_stanza_total( stanza_data = stanza_data,
                                  #Y_zz = Y_zz )
                                  Y_zz_g2 = Y_zz_g2 )
      Bnew_i[stanza_data$stanzainfo_s2z[,'s']] = Bnew_s2
    }

    # Average biomass
    for( i in seq_len(n_species) ){
      Bmean_ti[t,i] = mean(proj$y[,i])
    }

    # Record variables
    if( control$process_error == "epsilon" ){
      B_ti[t,] = Bnew_i
    }else{
      Bhat_ti[t,] = Bnew_i
      for( i in seq_len(n_species) ){
        if( !is.na(p$logtau_i[i]) ){
          B_ti[t,i] = out_initial$B_i[i] * exp(p$alpha_ti[t,i])
          epsilon_ti[t,i] = log( B_ti[t,i] / Bhat_ti[t,i] )
        }else{
          B_ti[t,i] = Bhat_ti[t,i]
          epsilon_ti[t,i] = 0
        }
      }
    }

    # Record other variables
    if( control$F_type=="integrated" ){
      Chat_ti[t,] = proj$y[nrow(proj$y),n_species+seq_len(n_species)]
    }else{
      Chat_ti[t,] = Bmean_ti[t,] * (1 - exp( -1 * exp(p$logF_ti[t,]) ))
    }
    DC_ij = as.matrix(p$DC_ij)
    M2_ti[t,] = (DC_ij %*% (Bmean_ti[t,] * exp(p$logQB_i))) / Bmean_ti[t,]

    # Record more using midpoint biomass Bmean_ti
    out_mean = dBdt( Time = 0,
                State = c(Bmean_ti[t,], rep(0,n_species)),
                Pars = p_t,
                what = "stuff" )

    # Must calculate during loop because G_ti is NA for t=1
    #tmp = adsparse_to_matrix(out$Q_ij)
    G_ti[t,] = out_mean$G_i
    g_ti[t,] = out_mean$g_i
    M_ti[t,] = out_mean$M_i
    m_ti[t,] = out_mean$m_i
    M2_ti[t,] = out_mean$M2_i
    m2_ti[t,] = out_mean$m2_i
    dBdt0_ti[t,] = out_mean$dBdt0_i
    Q_tij[t,,] = out_mean$Q_ij
    TL_ti[t,] = compute_tracer( Q_ij = out_mean$Q_ij,
                                inverse_method = control$inverse_method,
                                type_i = type_i,
                                tracer_i = rep(1,n_species) )

    # Compute trophic level at the end of each time
    # Not needed ... definiting stuff at annual averages
    #out = dBdt( Time = 0,
    #            State = c(B_ti[t,], rep(0,n_species)),
    #            Pars = p_t,
    #            what = "stuff" )
  }
  F_ti = exp(p$logF_ti)
  Z_ti = F_ti + M_ti 
  
  # likelihood
  Bobs_ti = OBS(Bobs_ti)
  Cobs_ti = OBS(Cobs_ti)
  Bexp_ti = B_ti * (rep(1,nrow(B_ti)) %*% t(exp(p$logq_i)))
  for( i in seq_len(n_species) ){
  for( t in seq_len(nrow(Bexp_ti)) ){
    if( !is.na(Bobs_ti[t,i]) ){
      loglik1_ti[t,i] = dnorm( log(Bobs_ti[t,i]), log(Bexp_ti[t,i]), exp(p$ln_sdB), log=TRUE)
    }
    if( (taxa %in% fit_eps)[i] ){
      loglik2_ti[t,i] = dnorm( epsilon_ti[t,i], 0, exp(p$logtau_i[i]), log=TRUE)
    }
    if( !is.na(Cobs_ti[t,i]) ){
      loglik3_ti[t,i] = dnorm( log(Cobs_ti[t,i]), log(Chat_ti[t,i]), exp(p$ln_sdC), log=TRUE)
    }
    if( (taxa %in% fit_nu)[i] ){
      loglik4_ti[t,i] = dnorm( p$nu_ti[t,i], 0, exp(p$logsigma_i[i]), log=TRUE)
    }
  }}
  #if(isFALSE(inherits(Bexp_ti,"advector"))) stop("Bexp_ti")

  # Comps
  dmultinomial = function( x, prob, log=TRUE ){
    r = lgamma(sum(x) + 1) + sum(x * log(prob+1e-20) - lgamma(x + 1))
    if(log){r}else{exp(r)}
  }
  # From gtools::ddirichlet
  ddirichlet <- function(x, alpha, log=TRUE) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s = (alpha - 1) * log(x)
    r = sum(s) - logD
    if(log){r}else{exp(r)}
  }
  #selex_index = 0
  for( index in seq_along(Nobs_ta_g2) ){
    g2 = match( names(Nobs_ta_g2)[index], settings$unique_stanza_groups )
    #which_z = which( stanza_data$X_zz[,'g2'] == g2 )
    #selex_a = plogis( (stanza_data$X_zz[which_z,'AGE'] - p$s50_z[index])/p$srate_z[index] )
    #which_a = which( stanza_data$X_zz[which_z,'AGE'] %in% unique(stanza_data$X_zz[which_z,'age_class']) )
    #selex_a = plogis( (stanza_data$X_zz[which_z[which_a],'AGE'] - p$s50_z[index])/p$srate_z[index] )
    Xg2_zz = stanza_data$X_zz_g2[[g2]]
    Yg2_tzz = Y_tzz_g2[[g2]]
    selex_a = plogis( (Xg2_zz[,'AGE'] - p$s50_z[index])/p$srate_z[index] )
    for( index2 in seq_len(nrow(Nobs_ta_g2[[index]])) ){
      t = match( rownames(Nobs_ta_g2[[index]])[index2], years )
      # Comps are average-year abundance (smears cohorts across adjacent years)
      #Nexp_a = rep(0, max(stanza_data$X_zz[which_z,'age_class']+1)) * p$s50_z / p$s50_z # 0 through MaxAge so +1 length
      #for(z in which_z){
      #  Nexp_a[stanza_data$X_zz[z,'age_class']+1] = Nexp_a[stanza_data$X_zz[z,'age_class']+1] + selex_a[z]*Y_tzz[t,z,'NageS']
      #}
      Nexp_a = rep(0, max(Xg2_zz[,'age_class']+1)) * p$s50_z[index] / p$s50_z[index] # 0 through MaxAge so +1 length
      for( z in seq_len(nrow(Xg2_zz)) ){
        Nexp_a[Xg2_zz[z,'age_class']+1] = Nexp_a[Xg2_zz[z,'age_class']+1] + selex_a[z] * exp(Yg2_tzz[t,z,'log_NageS'])
      }
      # Comps are end-of-year abundance
      #which_a = which( stanza_data$X_zz[which_z,'AGE'] %in% unique(stanza_data$X_zz[which_z,'age_class']) )
      #Nexp_a = selex_a * Y_tzz[t,which_z[which_a],'NageS']
      # Record comps
      Nexp_ta_g2[[index]][index2,] = Nexp_a[-1] + settings$min_agecomp_prob  # Remove age-0 ... add 1e-12 to avoid prob=0, which crashes gradients
      # Remove any NAs
      which_obs = which(!is.na((Nobs_ta_g2[[index]])[index2,]))
      obs = (Nobs_ta_g2[[index]])[index2,which_obs]
      prob = Nexp_ta_g2[[index]][index2,which_obs] / sum(Nexp_ta_g2[[index]][index2,which_obs],na.rm=TRUE)
      if(!any(is.na(obs))){
        if( settings$comp_weight == "multinom" ){
          loglik5_tg2[t,g2] = dmultinomial( obs, prob=prob, log=TRUE )
          # relative deviance: https://stats.stackexchange.com/questions/186560/what-is-multinomial-deviance-in-the-glmnet-package
          # dev5_tg2[t,g2] = loglik5_tg2[t,g2] - dmultinomial( obs, prob=obs/sum(obs), log=TRUE )
        }else if( settings$comp_weight == "dir" ){
          loglik5_tg2[t,g2] = ddirichlet( obs/sum(obs), alpha=prob*exp(p$compweight_z[index]), log=TRUE )
        }else{
          loglik5_tg2[t,g2] = ddirmult( obs, prob=prob, ln_theta=p$compweight_z[index], log=TRUE )
        }
      }
    }
  }

  # Empirical weight-at-age
  for( index in seq_along(Wobs_ta_g2) ){
    g2 = match( names(Wobs_ta_g2)[index], settings$unique_stanza_groups )
    #which_z = which( stanza_data$X_zz[,'g2'] == g2 )
    Xg2_zz = stanza_data$X_zz_g2[[g2]]
    Yg2_tzz = Y_tzz_g2[[g2]]
    #weight_index = max(weight_index) + 1:2
    for( index2 in seq_len(nrow(Wobs_ta_g2[[index]])) ){
      t = match( rownames(Wobs_ta_g2[[index]])[index2], years )
      #Wexp_a = Nexp_a = rep(0,max(stanza_data$X_zz[which_z,'age_class']+1)) # 0 through MaxAge so +1 length
      #for(z in which_z){
      #  Nexp_a[stanza_data$X_zz[z,'age_class']+1] = Nexp_a[stanza_data$X_zz[z,'age_class']+1] + Y_tzz[t,z,'NageS']
      #}
      #for(z in which_z){
      #  prop = Y_tzz[t,z,'NageS'] / Nexp_a[stanza_data$X_zz[z,'age_class']+1]
      #  Wexp_a[stanza_data$X_zz[z,'age_class']+1] = Wexp_a[stanza_data$X_zz[z,'age_class']+1] + prop * Y_tzz[t,z,'WageS']
      #}
      Wexp_a = Nexp_a = rep(0,max(Xg2_zz[,'age_class']+1)) # 0 through MaxAge so +1 length
      for( z in seq_len(nrow(Xg2_zz)) ){
        #Nexp_a[Xg2_zz[z,'age_class']+1] = Nexp_a[Xg2_zz[z,'age_class']+1] + Yg2_tzz[t,z,'NageS']
        Nexp_a[Xg2_zz[z,'age_class']+1] = Nexp_a[Xg2_zz[z,'age_class']+1] + exp(Yg2_tzz[t,z,'log_NageS'])
      }
      for( z in seq_len(nrow(Xg2_zz)) ){
        #prop = Yg2_tzz[t,z,'NageS'] / Nexp_a[Xg2_zz[z,'age_class']+1]
        prop = exp(Yg2_tzz[t,z,'log_NageS']) / Nexp_a[Xg2_zz[z,'age_class']+1]
        Wexp_a[Xg2_zz[z,'age_class']+1] = Wexp_a[Xg2_zz[z,'age_class']+1] + prop * Yg2_tzz[t,z,'WageS']
      }
      Wexp_ta_g2[[index]][index2,] = Wexp_a[-1] * exp(p$log_winf_z[index])   # Remove age-0
      #Wexp_ta_g2[[index]][index2,] = Wexp_a[-1]   # Remove age-0
      obs = (Wobs_ta_g2[[index]])[index2,]
      mu = (Wexp_ta_g2[[index]])[index2,]
      for( index3 in seq_along(obs) ){
        if( !is.na(obs[index3]) ){
          loglik6_tg2[t,g2] = loglik6_tg2[t,g2] + dnorm(log(obs[index3]), mean=log(mu[index3]), sd=exp(p$ln_sdW_z[index]), log=TRUE)
        }
      }
    }
  }

  for( g2 in seq_len(settings$n_g2) ){
  for( t in seq_len(nrow(Bexp_ti)) ){
    if( (settings$unique_stanza_groups %in% fit_phi)[g2] ){
      loglik7_tg2[t,g2] = dnorm( p$phi_tg2[t,g2], 0, exp(p$logpsi_g2[g2]), log=TRUE)
    }
  }}

  # Calculate priors
  log_prior_value = log_prior( p )

  # Remove NAs to deal with missing values in Bobs_ti and Cobs_ti
  jnll = jnll - ( sum(loglik1_ti) + sum(loglik2_ti) + sum(loglik3_ti) + sum(loglik4_ti) + sum(loglik5_tg2,na.rm=TRUE) + sum(loglik6_tg2) + sum(loglik7_tg2) + sum(log_prior_value,na.rm=TRUE) )
  
  ###############
  # Derived
  ###############

  # Steepness
  h_g2 = 0.2 * p$SpawnX_g2 / (p$SpawnX_g2 - 1 + 0.2)
  REPORT( h_g2 )

  # Abundance-at-age
  W_ta_g2 = N_ta_g2 = vector("list", length=length(Y_tzz_g2) )
  names(W_ta_g2) = names(N_ta_g2) = names(Y_tzz_g2 )
  for( g2 in seq_along(Y_tzz_g2) ){
    # Abundance-at-age
    N_ta2 = exp(Y_tzz_g2[[g2]][,,'log_NageS'])
    a_a2 = stanza_data$X_zz_g2[[g2]][,'age_class'] + 1
    N_at = apply( N_ta2, MARGIN=1, FUN=\(x) tapply(X=x, INDEX=a_a2, FUN=sum) )
    rownames(N_at) = unique(stanza_data$X_zz_g2[[g2]][,'age_class'])
    N_ta_g2[[g2]] = t(N_at)
    # Weight at age
    #W_at = N_at
    #W_ta_g2[[g2]] = t(N_at)
    for( t in seq_len(nrow(N_ta_g2[[g2]])) ){
      W_ta2 = Y_tzz_g2[[g2]][,,'WageS']
      N_ta2 = N_ta_g2[[g2]][,a_a2]
      prop_ta2 = exp(Y_tzz_g2[[g2]][,,'log_NageS']) / N_ta2
      #prop_a2 = exp(Y_tzz_g2[[g2]][t,,'log_NageS']) /
      W_ta2 = W_ta2 * prop_ta2
      W_at = apply( W_ta2, MARGIN=1, FUN=\(x) tapply(X=x, INDEX=a_a2, FUN=sum) )
      #W_at[,t] = tapply( W_a2, INDEX=a_a2, FUN=sum )
    }
    W_ta_g2[[g2]] = t(W_at)
  }
  REPORT( N_ta_g2 )
  REPORT( W_ta_g2 )

  ## weight-at-age for multistanza groups
  #for( g2 in seq_along(Y_tzz_g2) ){      # W_ta_g2[[g2]] =
  #  #N_ta_g2[[g2]] = matrix( 0, nrow=dim(Y_tzz_g2[[g2]])[1], ncol=length(unique(stanza_data$X_zz_g2[[g2]][,'age_class'])) )
  #  for( t in seq_len(nrow(N_ta_g2[[g2]])) ){
  #    # Abundance-at-age
  #    N_a2 = exp(Y_tzz_g2[[g2]][t,,'log_NageS'])
  #    a_a2 = stanza_data$X_zz_g2[[g2]][,'age_class'] + 1
  #    N_a = tapply( N_a2, INDEX=a_a2, FUN=sum )
  #    #N_ta_g2[[g2]][t,] = N_a
  #    # Abundance-weighted average weight-at-age
  #    W_a2 = Y_tzz_g2[[g2]][t,,'WageS']
  #    prop_a2 = N_a2 / N_a[a_a2]
  #    W_ta_g2[[g2]][t,] = tapply( W_a2 * prop_a2, INDEX=a_a2, FUN=sum )
  #  }
  #  colnames(N_ta_g2[[g2]]) = colnames(W_ta_g2[[g2]]) = unique(stanza_data$X_zz_g2[[g2]][,'age_class'])
  #}
  #REPORT( W_ta_g2 )

  # Allow for ADREPORT
  Wmat_g2 = p$Wmat_g2
  X_ij = 1 + exp(Xprime_ij)
  REPORT( Wmat_g2 )
  REPORT( X_ij )

  # Initial conditions
  REPORT( B0_i )
  REPORT( TL0_i )
  REPORT( G0_ti )
  REPORT( g0_ti )
  REPORT( M0_ti )
  REPORT( m0_ti )
  REPORT( M20_ti )
  REPORT( m20_ti )

  # Relative biomass
  BoverB0_ti = B_ti
  for(t in 1:nrow(BoverB0_ti)) BoverB0_ti[t,] = B_ti[t,] / B0_i
  REPORT( BoverB0_ti )

  # Reporting
  REPORT( B_ti )
  REPORT( TotalEggs_tg2 )
  REPORT( TotalSB_tg2 )
  REPORT( TotalZ_ts2 )
  if(control$process_error=="alpha") REPORT( Bhat_ti )
  REPORT( Chat_ti )
  REPORT( Bmean_ti )
  REPORT( out_initial )
  REPORT( Bexp_ti )
  REPORT( G_ti )
  REPORT( g_ti )
  REPORT( M_ti )
  REPORT( m_ti )
  REPORT( M2_ti )
  REPORT( m2_ti )
  REPORT( F_ti )
  REPORT( Z_ti )
  REPORT( Q_tij )
  REPORT( dBdt0_ti )
  REPORT( loglik1_ti )
  REPORT( loglik2_ti )
  REPORT( loglik3_ti )
  REPORT( loglik4_ti )
  REPORT( loglik5_tg2 )
  REPORT( loglik6_tg2 )
  REPORT( loglik7_tg2 )
  REPORT( log_prior_value )
  REPORT( jnll )
  REPORT( TL_ti )
  REPORT( Y_tzz_g2 )
  REPORT( stanza_data )
  REPORT( Nexp_ta_g2 )
  REPORT( Wexp_ta_g2 )

  if( settings$n_g2 >0 ){
    R0_g2 = p$baseR0_g2
    SB0_g2 = p$baseSB_g2
    REPORT( R0_g2 )
    REPORT( SB0_g2 )
  }

  if( length(control$derived_quantities) > 0 ){
    derived_values = mget( control$derived_quantities, ifnotfound=NA )
    REPORT( derived_values )
    ADREPORT( do.call("c", derived_values) )
  }

  return(jnll)
}
