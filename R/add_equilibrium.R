#' @title 
#' Compute equilibrium values
#'
#' @description 
#' Compute equilibrium values
#'
#' @param p list of parameters
#'
#' @details
#' todo
#'
#' @export
add_equilibrium <-
function( ecoparams,
          scale_solver,
          noB_i ) { 
  
  # Guidelines
  # no ifelse() or which() for advectors
  # rowSums(C_ij, na.rm=TRUE) seems to break RTMB
  # Don't use n_species, or redefine it locally (using former for now)
  
  # 
  B_i = exp(ecoparams$logB_i)
  QB_i = exp(ecoparams$logQB_i)
  PB_i = exp(ecoparams$logPB_i)
  EE_i = ecoparams$EE_i
  DC_ij = ecoparams$DC_ij
  
  # Get equilibrium values
  # Eq-1 of Lucey-etal-2020:  B_i*PB_i*EE_i = DC_ij%*%(B_i*QB_i)
  if( scale_solver == "joint" ){
    # Hardwire is.na = 0 as experiment ... doesn't cause advector lost attribute
    #noB_i = ifelse( is.na(B_i), 1, 0 )
    #noB_i = rep(0,n_species)
  
    # Use Rpath logic ... see Rpath_logic.R
    #browser()
    BioQB = B_i * QB_i
    C_ij  = DC_ij * ( rep(1,length(B_i)) %*% t(BioQB) ) # BioQB[col(DC_ij)] # ( rep(1,n_species) %*% t(BioQB) ) # 
    b_i = rowSums(C_ij[,which(noB_i==0)])    # NAs for missing B_i
    #b_i = ifelse( b_i==0, 1e-10, b_i )
    #b_i = c(17.5, 10.5, 2, 2, 0)
    
    # 
    #diag.a = ifelse( is.na(EE_i), B_i*PB_i, EE_i*PB_i )
    diag.a = B_i * PB_i
    diag.a[which(noB_i==1)] = EE_i[which(noB_i==1)] * PB_i[which(noB_i==1)]
    A = diag(diag.a)
    #A = matrix( c(90,0,0,0,0, 0,10.5,0,0,0, 0,0,0.2,0,0, 0,0,0,3,0, 0,0,0,0,0.1), nrow=5, byrow=TRUE)
    
    #
    QBDC = DC_ij * ( rep(1,length(B_i)) %*% t(QB_i) ) # QB_i[col(DC_ij)] # ( rep(1,n_species) %*% t(QB_i) ) # 
    # QBDC[is.na(QBDC)] = 0    # Not necessary given that all(!is.na(DC_ij))
    QBDCa = QBDC * ( rep(1,length(B_i)) %*% t(noB_i) ) # noB_i[col(as.matrix(QBDC))] # ( rep(1,n_species) %*% t(noB_i) ) # 
    A = A - QBDCa 
    
    # Generalized inverse does the actual solving
    #Invert A and multiple by b to get x (unknowns)
    # A is singular when some B_i are NA
    #x_i = MASS::ginv(A) %*% b_i
    #x_i = pseudoinverse(A) %*% b_i
    #x_i = solve(A, b_i) # solve(A) %*% b_i      
    x_i = solve(A) %*% b_i # solve(A) %*% b_i      
    #x_i = c(0.194, 1, 10, 0.6667, 0)
    
    # Substitute into vectors
    EE_i[which(noB_i==0)] = x_i[which(noB_i==0)]
    B_i[which(noB_i==1)] = x_i[which(noB_i==1)]
    # Using manual indices doesn't seem to fix things
    #EE_i[c(1,3,4,5)] = x_i[c(1,3,4,5)]
    #B_i[2] = x_i[2]
    #EE_i = c(0.194, 2.1, 10, 0.6667, 0 )
    #B_i = c(1, 1, 1, 1, 1)
  }else{
    # Derive EE_i from B_i 
    EE_i = (DC_ij %*% (B_i*QB_i))[,1] / ( PB_i * B_i )
  }
  # Equilibrium consumption
  Qe_ij = DC_ij * ( rep(1,length(B_i)) %*% t(B_i * QB_i) )
  # Growth efficiency, Lucy-2020 Eq. 2
  GE_i = PB_i / QB_i
  # Natural mortality
  M0_i = PB_i * (1 - EE_i)
  
  #
  new_list = list(
    Qe_ij = Qe_ij,
    GE_i = GE_i,
    M0_i = M0_i
  )
  ecoparams = c( ecoparams, new_list )
  ecoparams$logB_i = log(B_i)
  ecoparams$EE_i = EE_i
  return( ecoparams )
}
