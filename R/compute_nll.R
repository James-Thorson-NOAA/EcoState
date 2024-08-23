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
                       scale_solver = scale_solver,
                       noB_i = noB_i,
                       type_i = type_i )
  
  # Extract epsilon_ti (local copy be modified later)
  epsilon_ti = p$epsilon_ti
  if(process_error=="alpha"){
    epsilon_ti = matrix( 0, ncol=n_species, nrow=nrow(Bobs_ti) )
  }
  
  # unfished M0 and M2, and B_i solved for EE_i
  p_t = p
    p_t$epsilon_i = rep(0,n_species)
    p_t$logF_i = rep(-Inf,n_species)
  out_initial = dBdt( Time = 1, 
              State = c(p$logB_i,rep(0,n_species)), 
              Pars = p_t,
              what = "stuff" )
  
  # Objects to save
  TL_ti = dBdt0_ti = M_ti = m_ti = G_ti = g_ti = M2_ti = m2_ti = Bmean_ti = Chat_ti = B_ti = Bhat_ti = matrix( NA, ncol=n_species, nrow=nrow(Bobs_ti) )
  loglik1_ti = loglik2_ti = loglik3_ti = matrix( 0, ncol=n_species, nrow=nrow(Bobs_ti) )  # Missing = 0
  loglik4_tg2 = matrix( 0, nrow=nrow(Bobs_ti), ncol=length(settings$unique_stanza_groups) )
  Q_tij = array( NA, dim=c(nrow(Bobs_ti),n_species,n_species) )
  Nexp_ta_g2 = Nobs_ta_g2

  # Initial condition
  B_ti[1,] = out_initial$B_i * exp(p$delta_i)
  jnll = 0
  if( process_error == "alpha" ){
    epsilon_ti[1,] = p$alpha_ti[1,] 
  }
  Y_zz = p$Y_zz

  # 
  Y_tzz = array( 0.0, dim=c(nrow(Bobs_ti),dim(p$Y_zz)),
                 dimnames=list(NULL,rownames(p$Y_zz),colnames(p$Y_zz)) ) # dimnames-list dimension matching
  Y_tzz[1,,] = Y_zz

  # Loop through years
  for( t in 2:nrow(B_ti) ){
    # Assemble inputs
    p_t = p
    p_t$Y_zz = Y_zz
    p_t$logF_i = p$logF_ti[t,]

    # State-space or continuous innovations
    if( process_error == "epsilon" ){
      p_t$epsilon_i = epsilon_ti[t,]
    }else{
      p_t$epsilon_i = rep(0,n_species)
    }

    # RTMBode::ode requires y0 have names
    y0 = c(B_ti[t-1,], rep(0,n_species))
    names(y0) = paste0("var_",seq_along(y0))

    # Project dynamics
    #browser()
    proj = project_vars(
          f = dBdt,
          a = 0, 
          b = 1,
          n = n_steps,
          Pars = p_t,
          y0 = y0 )

    # Projec stanzas
    proj_stanzas = project_stanzas(
                 p = p_t,
                 stanza_data = stanza_data,
                 y = proj$y,
                 STEPS_PER_YEAR = settings$STEPS_PER_YEAR,
                 record_steps = FALSE,
                 correct_errors = TRUE )
    Y_zz = proj_stanzas$Y_zz
    Y_tzz[t,,] = Y_zz

    #
    Bnew_s2 = get_stanza_total( stanza_data = stanza_data,
                                Y_zz = p_t$Y_zz )
    Bnew_i = proj$y[nrow(proj$y),seq_len(n_species)]
    Bnew_i[p_t$stanzainfo_s2z[,'s']] = Bnew_s2

    # Average biomass
    for( i in seq_len(n_species) ){
      Bmean_ti[t,i] = mean(proj$y[,i])
    }

    # Record variables
    if( process_error == "epsilon" ){
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
    if( F_type=="integrated" ){
      Chat_ti[t,] = proj$y[nrow(proj$y),n_species+seq_len(n_species)]
    }else{
      Chat_ti[t,] = Bmean_ti[t,] * (1 - exp( -1 * exp(p$logF_ti[t,]) ))
    }
    DC_ij = as.matrix(DC_ij)
    M2_ti[t,] = (DC_ij %*% (Bmean_ti[t,] * exp(logQB_i))) / Bmean_ti[t,]

    # Record more using midpoint biomass Bmean_ti
    out = dBdt( Time = 0, 
                State = c(Bmean_ti[t,], rep(0,n_species)),
                Pars = p_t,
                what = "stuff" )

    # Must calculate during loop because G_ti is NA for t=1
    #tmp = adsparse_to_matrix(out$Q_ij)
    tmp = out$Q_ij
    G_ti[t,] = out$G_i
    g_ti[t,] = out$g_i
    M_ti[t,] = out$M_i
    m_ti[t,] = out$m_i
    M2_ti[t,] = out$M2_i
    m2_ti[t,] = out$m2_i
    Q_tij[t,,] = tmp
    dBdt0_ti[t,] = out$dBdt0_i
    # Compute trophic level
    TL_ti[t,] = compute_tracer( Q_ij = tmp,
                                inverse_method = inverse_method,
                                type_i = type_i,
                                tracer_i = rep(1,n_species) )
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
    if( !is.na(p$logtau_i[i]) ){
      loglik2_ti[t,i] = dnorm( epsilon_ti[t,i], 0, exp(p$logtau_i[i]), log=TRUE)
    }
    if( !is.na(Cobs_ti[t,i]) ){
      loglik3_ti[t,i] = dnorm( log(Cobs_ti[t,i]), log(Chat_ti[t,i]), exp(p$ln_sdC), log=TRUE)
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
  selex_index = 0
  for( index in seq_along(Nobs_ta_g2) ){
    g2 = match( names(Nobs_ta_g2)[index], settings$unique_stanza_groups )
    which_z = which( stanza_data$X_zz[,'g2'] == g2 )
    selex_index = max(selex_index) + seq_len( switch(settings$comp_weight,"multinom"=2,"dir"=3,"dirmult"=3) ) # CHANGE WITH NUMBER OF PARAMETERS
    selex_pars = p$selex_z[selex_index]
    selex_a = plogis( (stanza_data$X_zz[which_z,'AGE'] - selex_pars[1])/selex_pars[2] )
    for( index2 in seq_len(nrow(Nobs_ta_g2[[index]])) ){
      t = match( rownames(Nobs_ta_g2[[index]])[index2], years )
      Nexp_a = rep(0,max(stanza_data$X_zz[which_z,'age_class']+1)) # 0 through MaxAge so +1 length
      for(z in which_z){
        Nexp_a[stanza_data$X_zz[z,'age_class']+1] = Nexp_a[stanza_data$X_zz[z,'age_class']+1] + selex_a[z]*Y_tzz[t,z,'NageS']
      }
      Nexp_ta_g2[[index]][index2,] = Nexp_a[-1] / sum(Nexp_a[-1],na.rm=TRUE)  # Remove age-0
      obs = (Nobs_ta_g2[[index]])[index2,]
      prob = Nexp_ta_g2[[index]][index2,]
      if( settings$comp_weight == "multinom" ){
        loglik4_tg2[t,g2] = dmultinomial( obs, prob=prob, log=TRUE )
      }else if( settings$comp_weight == "dir" ){
        loglik4_tg2[t,g2] = ddirichlet( obs/sum((Nobs_ta_g2[[index]])[index2,]), alpha=prob * exp(selex_pars[3]), log=TRUE )
      }else{
        #browser()
        loglik4_tg2[t,g2] = ddirmult( obs, prob=prob, ln_theta=selex_pars[3], log=TRUE )
      }
    }
  }

  # Remove NAs to deal with missing values in Bobs_ti and Cobs_ti
  jnll = jnll - ( sum(loglik1_ti) + sum(loglik2_ti) + sum(loglik3_ti) + sum(loglik4_tg2) )
  
  # Reporting
  REPORT( B_ti )
  if(process_error=="alpha") REPORT( Bhat_ti )
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
  REPORT( loglik4_tg2 )
  REPORT( jnll )
  REPORT( TL_ti )
  REPORT( Y_tzz )
  REPORT( stanza_data )
  REPORT( Nexp_ta_g2 )

  if( sdreport_detail >= 1 ){
    ADREPORT( B_ti )
  }
  if( sdreport_detail >= 2 ){
    ADREPORT( Chat_ti )
    ADREPORT( Bexp_ti )
    ADREPORT( G_ti )
    ADREPORT( g_ti )
    ADREPORT( M_ti )
    ADREPORT( m_ti )
    ADREPORT( F_ti )
    ADREPORT( Z_ti )
    ADREPORT( dBdt0_ti )
    ADREPORT( TL_ti )
  }

  return(jnll)
}
