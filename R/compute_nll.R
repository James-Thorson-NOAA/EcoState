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
  
  # Objects to save
  dBdt0_ti = deltaBB_ti = P_ti = M_ti = G_ti = M2_ti = Bhatmean_ti = Chat_ti = Bhat_ti = matrix( NA, ncol=n_species, nrow=nrow(Bobs_ti) )
  Q_tij = array( NA, dim=c(nrow(Bobs_ti),n_species,n_species) )
  
  # Initial condition
  Bhat_ti[1,] = exp(p$logB_i + p$logB0ratio_i)
  
  if( FALSE ){
    t = 2
    p_t = p
    p_t$deltaB_i = p$deltaB_ti[t,]
    p_t$logF_i = p$logF_ti[t,]
    n_steps = 10
    y0 = c(Bhat_ti[t-1,], rep(0,n_species))
    proj = myode(
          f = dBdt,
          a = 0, 
          b = 1,
          n = n_steps,
          Pars = p_t,
          y0 = y0 )
  }
  
  # Loop through years
  for( t in 2:nrow(Bhat_ti) ){
    # Assemble inputs
    p_t = p
    p_t$deltaB_i = p$deltaB_ti[t,]
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
    Chat_ti[t,] = proj$y[nrow(proj$y),n_species+seq_len(n_species)]
    M2_ti[t,] = (DC_ij %*% (Bhatmean_ti[t,] * exp(logQB_i))) / Bhatmean_ti[t,]

    # Record more using midpoint biomass Bhatmean_ti
    out = dBdt( Time = 0, 
                State = c(Bhatmean_ti[t,], rep(0,n_species)),
                Pars = p_t,
                what = "stuff" )
    G_ti[t,] = out$G_i
    M_ti[t,] = out$M_i
    Q_tij[t,,] = out$Q_ij
    dBdt0_ti[t,] = out$dBdt0_i
    # Must calculate during loop because G_ti is NA for t=1
    P_ti[t,] = G_ti[t,] / Bhat_ti[t,]
    deltaBB_ti[t,] = p$deltaB_ti[t,]
  }
  F_ti = exp(p$logF_ti)
  Z_ti = F_ti + M_ti 
  
  # likelihood
  Bexp_ti = Bhat_ti * (rep(1,nrow(Bhat_ti)) %*% t(exp(p$logq_i)))
  jnll = 0
  for( i in seq_len(n_species) ){
    jnll = jnll - sum( dnorm( log(Bobs_ti[,i]), log(Bexp_ti[,i]), exp(p$ln_sdB), log=TRUE), na.rm=TRUE )
    if( !is.na(p$logtau_i[i]) ){
      jnll = jnll - sum( dnorm(p$deltaB_ti[,i], 0, exp(p$logtau_i[i]), log=TRUE) )
    }
    jnll = jnll - sum( dnorm(log(Cobs_ti[,i]), log(Chat_ti[,i]), exp(p$ln_sdC), log=TRUE), na.rm=TRUE )
  }
  
  # unfished M0 and M2
  out = dBdt( Time = 1, 
        State = c(p$logB_i,rep(0,n_species)), 
        Pars = p_t,
        what = "stuff" )
  
  # Reporting
  REPORT( Bhat_ti )
  REPORT( Chat_ti )
  REPORT( M2_ti )
  REPORT( Bhatmean_ti )
  REPORT( out )
  REPORT( Bexp_ti )
  REPORT( G_ti )
  REPORT( M_ti )
  REPORT( F_ti )
  REPORT( Z_ti )
  REPORT( P_ti )
  REPORT( deltaBB_ti )
  REPORT( Q_tij )
  REPORT( dBdt0_ti );
  ADREPORT( Bhat_ti )
  ADREPORT( Chat_ti )
  ADREPORT( Bexp_ti )
  ADREPORT( G_ti )
  ADREPORT( M_ti )
  ADREPORT( F_ti )
  ADREPORT( Z_ti )
  ADREPORT( P_ti )
  ADREPORT( deltaBB_ti )
  ADREPORT( dBdt0_ti )
  
  return(jnll)
}
