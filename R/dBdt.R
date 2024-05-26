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
#' @export
dBdt <-
function( Time, 
          State, 
          Pars, 
          what="dBdt" ){     

  # Inputs
  getAll(Pars)
  Bt_i = State[1:n_species]
  B_i = exp(logB_i)
  QB_i = exp(logQB_i)
  PB_i = exp(logPB_i)
  V_ij = exp(logV_ij)

  # Get equilibrium values
  # Derive EE_i from B_i, Eq-1 of Lucey-etal-2020
  # B_i*PB_i*EE_i = DC_ij%*%(B_i*QB_i)
  EE_i = (DC_ij %*% (B_i*QB_i))[,1] / ( PB_i * B_i )
  # Equilibrium consumption
  Qe_ij = DC_ij * ( rep(1,n_species) %*% t(B_i * QB_i) )
  # Growth efficiency, Lucy-2020 Eq. 2
  GE_i = PB_i / QB_i
  # Natural mortality
  M0_i = PB_i * (1 - EE_i)

  # Compute dyynamics 
  # Predator and prey abundance relative to equilibrium
  Ypred_j = Bt_i / B_i
  Ypred_ij = rep(1,n_species) %*% t(Ypred_j)
  Yprey_i = Bt_i / B_i
  # Consumption = Equilibrium * Pred_functional_response * Prey_functional_response
  Q_ij = Qe_ij * ( V_ij * Ypred_ij / ( V_ij - 1 + Ypred_ij ) ) * Yprey_i
  # Calculate growth G_i (called C_i originally but conflicts with catch C_ti)
  G_i = GE_i * colSums(Q_ij)
  # Replace production for consumption for primary producers, including self-limitation via V_ij
  G_i = ifelse( colSums(DC_ij)==0, PB_i * Bt_i, G_i )
  #G_i = ifelse( colSums(DC_ij)==0, ( diag(V_ij) * Ypred_j / ( diag(V_ij) - 1 + Ypred_j ) ) * B_i * PB_i, G_i )
  
  # Assemble dynamics
  dBdt0_i = G_i - rowSums(Q_ij) - M0_i*Bt_i
  # Include stochasticity ... as function of Bt_i
  dBdt_i = dBdt0_i + deltaB_i*Bt_i
  # Augment with fishing mortality and catches
  dBdt_i = c( dBdt_i - Bt_i*exp(logF_i), Bt_i*exp(logF_i) )
  
  # Outputs
  if(what=="dBdt"){
    return(dBdt_i)
  }else{
    # Calculate mortality
    M_i = rowSums(Q_ij)/Bt_i + M0_i
    # Predation mortality .. removed because it doesn't vary over time
    #M2_i = (DC_ij %*% (B_i*QB_i))[,1] / B_i
    # Bundle and return
    return( list(EE_i=EE_i, GE_i=GE_i, M0_i=M0_i, Q_ij=Q_ij, G_i=G_i, M_i=M_i, Qe_ij=Qe_ij, dBdt0_i=dBdt0_i) )
  }
}
