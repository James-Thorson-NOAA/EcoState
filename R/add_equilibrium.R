#' @title 
#' Compute equilibrium values
#'
#' @inheritParams ecostate
#' @inheritParams ecostate_control
#'
#' @description 
#' Compute equilibrium values
#'
#' @param ecoparams list of parameters
#'
#' @details
#' todo
#'
#' @export
add_equilibrium <-
function( ecoparams,
          scale_solver,
          noB_i,
          type_i ) { 
  
  # Guidelines
  # no ifelse() or which() for advectors
  # rowSums(C_ij, na.rm=TRUE) seems to break RTMB
  # Don't use n_species, or redefine it locally (using former for now)
  
  # Indicators 
  which_primary = which( type_i=="auto" )
  which_detritus = which( type_i=="detritus" )
  
  # 
  B_i = exp(ecoparams$logB_i)
  QB_i = exp(ecoparams$logQB_i)
  PB_i = exp(ecoparams$logPB_i)
  EE_i = ecoparams$EE_i
  DC_ij = ecoparams$DC_ij
  U_i = ecoparams$U_i
  
  # Get equilibrium values
  # Eq-1 of Lucey-etal-2020:  B_i*PB_i*EE_i = DC_ij%*%(B_i*QB_i)
  if( scale_solver == "joint" ){
    # Hardwire is.na = 0 as experiment ... doesn't cause advector lost attribute
    #noB_i = ifelse( is.na(B_i), 1, 0 )
    #noB_i = rep(0,n_species)
  
    # Use Rpath logic ... see Rpath_logic.R
    #browser()
    QB_i[c(which_primary,which_detritus)] = 0
    BioQB = B_i * QB_i
    C_ij  = DC_ij * ( rep(1,length(B_i)) %*% t(BioQB) ) # BioQB[col(DC_ij)] # ( rep(1,n_species) %*% t(BioQB) ) # 
    b_i = rowSums(C_ij[,which(noB_i==0),drop=FALSE])    # NAs for missing B_i
    #b_i = ifelse( b_i==0, 1e-10, b_i )
    #b_i = c(17.5, 10.5, 2, 2, 0)
    
    # 
    diag.a = B_i * PB_i
    diag.a[which(noB_i==1)] = EE_i[which(noB_i==1)] * PB_i[which(noB_i==1)]
    # diag(vector) doesn't work when length(vector)=1, so using explicit construction
    A = matrix(0, nrow=length(B_i), ncol=length(B_i))
    A[cbind(seq_along(B_i),seq_along(B_i))] = diag.a
    
    #
    QBDCa = DC_ij * ( rep(1,length(B_i)) %*% t(QB_i*noB_i) ) # QB_i[col(DC_ij)] # ( rep(1,n_species) %*% t(QB_i) ) # 
    # QBDC[is.na(QBDC)] = 0    # Not necessary given that all(!is.na(DC_ij))
    A = A - QBDCa 
    
    # Generalized inverse does the actual solving
    #Invert A and multiple by b to get x (unknowns)
    # A is singular when some B_i are NA
    #x_i = MASS::ginv(A) %*% b_i
    #x_i = pseudoinverse(A) %*% b_i
    #x_i = solve(A, b_i) # solve(A) %*% b_i      
    x_i = solve(A) %*% b_i # solve(A) %*% b_i      
    
    # Substitute into vectors
    EE_i[which(noB_i==0)] = x_i[which(noB_i==0)]
    B_i[which(noB_i==1)] = x_i[which(noB_i==1)]
  }else{
    # Derive EE_i from B_i 
    QB_i[c(which_primary,which_detritus)] = 0
    EE_i = (DC_ij %*% (B_i*QB_i))[,1] / ( PB_i * B_i )
    # Sanity check
    # (PB_i * EE_i * B_i) = (DC_ij %*% (B_i*QB_i))[,1]
  }
  # Equilibrium consumption
  Qe_ij = DC_ij * ( rep(1,length(B_i)) %*% t(B_i * QB_i) )
  # Growth efficiency, Lucy-2020 Eq. 2
  GE_i = PB_i / QB_i
  # Natural mortality
  m0_i = PB_i * (1 - EE_i)
  
  # Calculate detritus turnover
  detritus_input = sum(Qe_ij %*% matrix(U_i)) + sum( m0_i * B_i )
  detritus_output = sum(Qe_ij[which_detritus,])
  detritus_turnover = (detritus_input-detritus_output) / B_i[which_detritus]        # detrius_input - detritus_output - B*detritus_turnover = 0
  
  #
  new_list = list(
    Qe_ij = Qe_ij,
    GE_i = GE_i,
    m0_i = m0_i,
    detritus_turnover = detritus_turnover
  )
  ecoparams = c( ecoparams, new_list )
  # Overwrite values
  ecoparams$logB_i = log(B_i)
  ecoparams$EE_i = EE_i
  return( ecoparams )
}
