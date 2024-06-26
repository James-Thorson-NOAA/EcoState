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
          what = "dBdt" ){     

  # Inputs from function call
  getAll(Pars)
  Bt_i = State[1:n_species]

  # Inputs from local environment
  QB_i = exp(logQB_i)
  PB_i = exp(logPB_i)
  V_ij = exp(Vprime_ij) + 1
  B_i = exp(logB_i)
  # EE_i is in Pars
  # U_i is in Pars
  F_i = exp(logF_i)
  
  # Compute dyynamics 
  # Predator and prey abundance relative to equilibrium
  Ypred_j = Bt_i / B_i
  Ypred_ij = rep(1,n_species) %*% t(Ypred_j)
  Yprey_i = Bt_i / B_i
  Yprey_ij = Yprey_i %*% t(rep(1,n_species))
  U_ij = rep(1,n_species) %*% t(U_i)
  # Consumption = Equilibrium * Pred_functional_response * Prey_functional_response
  Q_ij = Qe_ij * ( V_ij * Ypred_ij / ( V_ij - 1 + Ypred_ij ) ) * Yprey_ij

  # Calculate growth G_i (called C_i originally but conflicts with catch C_ti)
  G_i = GE_i * colSums(Q_ij)
  # Replace production for consumption for primary producers, including self-limitation via V_ij
  numerator = diag(V_ij[which_primary,which_primary,drop=FALSE]) * Yprey_i[which_primary]
  denominator = ( diag(V_ij[which_primary,which_primary,drop=FALSE]) - 1 + Yprey_i[which_primary] )
  G_i[which_primary] = PB_i[which_primary] * B_i[which_primary] * numerator / denominator
  # Replace for detritus
  G_i[which_detritus] = sum(Q_ij * U_ij) + sum(m0_i * Bt_i)
  
  # Mortalities
  M0_i = m0_i * Bt_i
  M2_i = rowSums(Q_ij)
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
    #M2_i = (DC_ij %*% (B_i*QB_i))[,1] / B_i
    # Bundle and return
    Return = list( B_i = B_i, 
                   EE_i = EE_i, 
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
