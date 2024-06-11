
#' @title Calculate trophic level
#'
#' @description Calculate trophic level from consumption matrix
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
#' @param inverse_method whether to use pseudoinverse or standard inverse
#'
#' @details
#' Trophic level \eqn{l_i} for each predator \eqn{i} is defined as:
#'
#' \deqn{ \mathbf{l - 1 = l Q^*} }
#'
#' where \eqn{\mathbf{Q*}} is the proportion consumption for each predator (column)
#' of different prey (rows).  We identify primary producers as any taxa with no
#' consumption (a column of 0s), and assign them as the first trophic level.
#'
#' @export
get_trophic_level <-
function( Q_ij,
          inverse_method = c("Penrose_moore", "Standard"),
          which_primary = which(colSums(Q_ij)==0) ){
  # Start
  inverse_method = match.arg(inverse_method)
  
  # Rescale consumption
  colsums = colSums(Q_ij)
  B = diag(1/colsums)
  #inverse_denom = rep(1,nrow(Q_ij)) %*% t(1/colsums)
  Qprime_ij = Q_ij %*% B

  # Identify primary producers
  Qprime_ij[,which_primary] = 0
  
  # Solve for trophic level
  if( inverse_method == "Penrose_moore" ){
    inverse_IminusQ = corpcor::pseudoinverse( diag(nrow(Qprime_ij)) - Qprime_ij )
  }else{
    inverse_IminusQ = solve( diag(nrow(Qprime_ij)) - Qprime_ij )
  }
  l_i = t(rep(1,nrow(Qprime_ij))) %*% inverse_IminusQ
  return( l_i )
}
