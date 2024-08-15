#' @title Print EcoSim parameters
#'
#' @export
convert_rpath <- 
function( rpath_params,
          which_keep ){

  # Exclude fisheries
  if(missing(which_keep)){
    which_keep = which( rpath_params$model$Type != 3 )
  }
  
  # Define
  taxa = rpath_params$model$Group[which_keep]
  
  #stgroups <- c( "JuvRoundfish1"="round1", "AduRoundfish1"="round1", 
  #               "JuvRoundfish2"="round2",  "AduRoundfish2"="round2",  
  #               "JuvFlatfish1"="flat1", "AduFlatfish1"="flat1",
  #               "JuvFlatfish2"="flat2", "AduFlatfish2"="flat2" )
  rpath_params$stanzas$stindiv
  rpath_params$stanzas$stgroups
  stgroups = unlist(rpath_params$stanzas$stgroups[,'StanzaGroup'])[unlist(rpath_params$stanzas$stindiv[,'StGroupNum'])]
  names(stgroups) = unlist(rpath_params$stanzas$stindiv[,'Group'])
  
  # Convert type
  type  <- sapply( rpath_params$model$Type[which_keep]+1, 
                   FUN = switch,
                   "hetero", "auto", "detritus" )
  
  # Extract stanza stuff
  K = unlist(rpath_params$stanzas$stgroups[,'VBGF_Ksp']) 
  Wmat = unlist(rpath_params$stanzas$stgroups[,'Wmat'])  
  d = unlist(rpath_params$stanzas$stgroups[,'VBGF_d'])
  
  # Convert Amax
  Amax = ceiling(unlist(rpath_params$stanzas$stindiv[,'Last'] + 1) / 12)
  
  # Extract Ecopath vectors
  PB = unlist(rpath_params$model[seq_along(taxa),'PB'])
  QB = unlist(rpath_params$model[seq_along(taxa),'QB'])
  QB = ifelse( is.na(QB), PB/unlist(rpath_params$model[seq_along(taxa),'ProdCons']), QB )
  B = unlist(rpath_params$model[seq_along(taxa),'Biomass'])
  EE = unlist(rpath_params$model[seq_along(taxa),'EE'])
  U = unlist(rpath_params$model[seq_along(taxa),'Unassim'])
  
  # Extract Ecopath matrices
  DC = cbind( as.matrix(rpath_params$diet[seq_along(taxa),-1]), Detritus = 0 )
  DC = ifelse( is.na(DC), 0, DC )
  X = array( 2, dim=c(length(taxa),length(taxa)) )
  SpawnX = rep( 2, length(unique(stgroups)) )
  
  # Label names
  names(K) = names(Wmat) = names(d) = names(SpawnX) = unique(stgroups)
  names(Amax) = names(stgroups)
  names(type) = names(PB) = names(QB) = names(B) = names(EE) = names(U) = taxa
  dimnames(DC) = dimnames(X) = list(taxa,taxa)

  # 
  out = list(
    taxa = taxa,
    stgroups = stgroups,
    type = type,
    PB = PB,
    QB = QB,
    B = B,
    EE = EE,
    U = U,
    DC = DC,
    X = X,
    SpawnX = SpawnX,
    K = K,
    Wmat = Wmat,
    d = d,
    Amax = Amax
  )
  return(out)
}
