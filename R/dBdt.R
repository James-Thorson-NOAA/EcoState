#' @title 
#' Dynamics from EcoSim
#'
#' @description 
#' Compute system of differential equations representing EcoSim dynamics
#'
#' @param Time todo
#' @param State todo
#' @param Pars todo
#' @param what todo
#'
#' @details
#' todo
#'
#' @return
#' An object (list) of ranges. Elements include:
#' \describe{
#' \item{G_i}{Biomass growth per time}
#' \item{g_i}{Biomass growth per time per biomass}
#' \item{M2_i}{Consumptive mortality per time}
#' \item{m2_i}{Consumptive mortality per time per biomass}
#' \item{M_i}{Natural mortality per time}
#' \item{m_i}{Natural mortality per time per biomass (i.e., m2_i + m0_i)}
#' \item{Q_ij}{Consumption per time for prey (rows) by predator (columns)}
#' }
#'
#' @export
dBdt <-
function( Time, 
          State, 
          Pars,
          F_type = "integrated",
          what = "dBdt" ){     

  # Inputs from function call
  RTMB::getAll(Pars)
  Bt_i = State[1:n_species]

  # Inputs from local environment
  QB_i = exp(logQB_i)
  PB_i = exp(logPB_i)
  X_ij = exp(Xprime_ij) + 1
  B_i = exp(logB_i)
  # EE_i is in Pars
  # U_i is in Pars
  F_i = exp(logF_i)
  
  # passed via environment:  type_i
  
  # Indicators
  which_primary = which( type_i=="auto" )
  which_detritus = which( type_i=="detritus" )
  
  # Compute dyynamics 
  # Predator and prey abundance relative to equilibrium
  Ypred_j = Bt_i / B_i
  Ypred_ij = rep(1,n_species) %*% t(Ypred_j)
  Yprey_i = Bt_i / B_i
  Yprey_ij = Yprey_i %*% t(rep(1,n_species))
  nu_ij = rep(1,n_species) %*% t(nu_i)
  # Consumption = Equilibrium * Pred_functional_response * Prey_functional_response
  #Qe_ij = adsparse_to_matrix(Qe_ij)
  Q_ij = Qe_ij * ( X_ij * Ypred_ij / ( X_ij - 1 + Ypred_ij ) ) * Yprey_ij * exp(nu_ij)
  #Q_ij = elementwise_product( Qe_ij, ( X_ij * Ypred_ij / ( X_ij - 1 + Ypred_ij ) ) * Yprey_ij )

  # Calculate growth G_i (called C_i originally but conflicts with catch C_ti)
  G_i = GE_i * colSums(Q_ij)
  #G_i = GE_i * (t(rep(1,n_species)) %*% Q_ij)[1,]
  # Replace production for consumption for primary producers, including self-limitation via X_ij
  numerator = X_ij[cbind(which_primary,which_primary)] * Yprey_i[which_primary]
  denominator = ( X_ij[cbind(which_primary,which_primary)] - 1 + Yprey_i[which_primary] )
  G_i[which_primary] = PB_i[which_primary] * B_i[which_primary] * numerator / denominator
  # Replace for detritus
  #G_i[which_detritus] = sum(Q_ij * U_ij) + sum(m0_i * Bt_i)
  G_i[which_detritus] = sum(Q_ij %*% matrix(U_i)) + sum(m0_i * Bt_i)
  
  # Mortalities
  M0_i = m0_i * Bt_i
  M2_i = rowSums(Q_ij)
  #M2_i = (Q_ij %*% rep(1,n_species))[,1]
  m2_i = M2_i / Bt_i
  
  # Total mortality
  m_i = m2_i + m0_i
  # Replace for detritus
  m_i[which_detritus] = m2_i[which_detritus] + detritus_turnover
  M_i = m_i * Bt_i

  # Assemble dynamics
  dBdt0_i = G_i - M_i
  # Include stochasticity ... as function of Bt_i
  dBdt1_i = dBdt0_i + epsilon_i*Bt_i
  # Augment with fishing mortality and catches
  if( F_type=="integrated" ){
    dBdt_i = c( dBdt1_i - Bt_i*F_i, Bt_i*F_i )
  }else{
    dBdt_i = dBdt1_i - Bt_i*F_i
  }

  # Outputs
  if(what=="dBdt"){
    Return = dBdt_i
  }else{
    # Convert to rates
    g_i = G_i / Bt_i

    # Predation mortality .. removed because it doesn't vary over time
    # Bundle and return
    Return = list( B_i = B_i, 
                   EE_i = EE_i,
                   PB_i = PB_i, 
                   QB_i = QB_i, # For filled in stanza values 
                   GE_i = GE_i, 
                   m0_i = m0_i, 
                   M0_i = M0_i, 
                   Q_ij = Q_ij, 
                   G_i = G_i,
                   g_i = g_i, 
                   M_i = M_i, 
                   m_i = m_i, 
                   M2_i = M2_i, 
                   m2_i = m2_i, 
                   Qe_ij = Qe_ij,
                   detritus_turnover = detritus_turnover, 
                   dBdt0_i = dBdt0_i, 
                   dBdt1_i = dBdt1_i, 
                   dBdt_i = dBdt_i)
  }
  return( Return )
}
