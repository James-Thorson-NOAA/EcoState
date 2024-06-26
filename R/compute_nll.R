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
  
  p = add_equilibrium( p,
                       scale_solver = scale_solver,
                       noB_i = noB_i )
  
  # unfished M0 and M2, and B_i solved for EE_i
  p_t = p
    p_t$epsilon_i = rep(0,n_species)
    p_t$logF_i = rep(-Inf,n_species)
  out_initial = dBdt( Time = 1, 
              State = c(p$logB_i,rep(0,n_species)), 
              Pars = p_t,
              what = "stuff" )
  
  # Objects to save
  TL_ti = dBdt0_ti = M_ti = m_ti = G_ti = g_ti = M2_ti = m2_ti = Bhatmean_ti = Chat_ti = Bhat_ti = matrix( NA, ncol=n_species, nrow=nrow(Bobs_ti) )
  loglik1_ti = loglik2_ti = loglik3_ti = matrix( 0, ncol=n_species, nrow=nrow(Bobs_ti) )  # Missing = 0
  Q_tij = array( NA, dim=c(nrow(Bobs_ti),n_species,n_species) )
  
  # Initial condition
  Bhat_ti[1,] = out_initial$B_i * exp(p$delta_i)
  jnll = 0

  # Loop through years
  for( t in 2:nrow(Bhat_ti) ){
    # Assemble inputs
    p_t = p
    p_t$epsilon_i = p$epsilon_ti[t,]
    p_t$logF_i = p$logF_ti[t,]

    # RTMBode::ode requires y0 have names
    y0 = c(Bhat_ti[t-1,], rep(0,n_species))
    names(y0) = paste0("var_",seq_along(y0))

    # Project dynamics
    proj = project_vars(
          f = dBdt,
          a = 0, 
          b = 1,
          n = n_steps,
          Pars = p_t,
          y0 = y0 )

    # Average biomass
    for( i in seq_len(n_species) ){
      Bhatmean_ti[t,i] = mean(proj$y[,i])
    }

    # Record variables
    Bhat_ti[t,] = proj$y[nrow(proj$y),seq_len(n_species)]
    if( F_type=="integrated" ){
      Chat_ti[t,] = proj$y[nrow(proj$y),n_species+seq_len(n_species)]
    }else{
      Chat_ti[t,] = Bhatmean_ti[t,] * (1 - exp( -1 * exp(p$logF_ti[t,]) ))
    }
    M2_ti[t,] = (DC_ij %*% (Bhatmean_ti[t,] * exp(logQB_i))) / Bhatmean_ti[t,]

    # Record more using midpoint biomass Bhatmean_ti
    out = dBdt( Time = 0, 
                State = c(Bhatmean_ti[t,], rep(0,n_species)),
                Pars = p_t,
                what = "stuff" )
    # Must calculate during loop because G_ti is NA for t=1
    G_ti[t,] = out$G_i
    g_ti[t,] = out$g_i
    M_ti[t,] = out$M_i
    m_ti[t,] = out$m_i
    M2_ti[t,] = out$M2_i
    m2_ti[t,] = out$m2_i
    Q_tij[t,,] = out$Q_ij
    dBdt0_ti[t,] = out$dBdt0_i
    # Compute trophic level
    TL_ti[t,] = compute_tracer( out$Q_ij, 
                                inverse_method = inverse_method,
                                which_primary = which_primary,
                                tracer_i = rep(1,n_species) )
  }
  F_ti = exp(p$logF_ti)
  Z_ti = F_ti + M_ti 
  
  # likelihood
  Bexp_ti = Bhat_ti * (rep(1,nrow(Bhat_ti)) %*% t(exp(p$logq_i)))
  for( i in seq_len(n_species) ){
  for( t in seq_len(nrow(Bexp_ti)) ){
    if( !is.na(Bobs_ti[t,i]) ){
      loglik1_ti[t,i] = dnorm( log(Bobs_ti[t,i]), log(Bexp_ti[t,i]), exp(p$ln_sdB), log=TRUE)
    }
    if( !is.na(p$logtau_i[i]) ){
      loglik2_ti[t,i] = dnorm( p$epsilon_ti[t,i], 0, exp(p$logtau_i[i]), log=TRUE)
    }
    if( !is.na(Cobs_ti[t,i]) ){
      loglik3_ti[t,i] = dnorm( log(Cobs_ti[t,i]), log(Chat_ti[t,i]), exp(p$ln_sdC), log=TRUE)
    }
  }}
  # Remove NAs to deal with missing values in Bobs_ti and Cobs_ti
  jnll = jnll - ( sum(loglik1_ti) + sum(loglik2_ti) + sum(loglik3_ti) ) 
  
  # Reporting
  REPORT( Bhat_ti )
  REPORT( Chat_ti )
  REPORT( Bhatmean_ti )
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
  REPORT( jnll )
  REPORT( TL_ti )

  ADREPORT( Bhat_ti )
  ADREPORT( Chat_ti )
  ADREPORT( Bexp_ti )
  ADREPORT( G_ti )
  ADREPORT( g_ti )
  ADREPORT( M_ti )
  ADREPORT( m_ti )
  ADREPORT( F_ti )
  ADREPORT( Z_ti )
  ADREPORT( dBdt0_ti )
  #if( inverse_method=="Standard" ){
    ADREPORT( TL_ti )
  #}
  
  return(jnll)
}
