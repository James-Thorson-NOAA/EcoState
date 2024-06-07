
#' @title Calculate trophic level
#'
#' @description Calculate trophic level from consumption matrix
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
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
get_trophic_level = function( Q_ij ){
  # Rescale consumption
  Qprime_ij = Q_ij / outer(rep(1,nrow(Q_ij)),colSums(Q_ij))

  # Identify primary producers
  which_pp = which(colSums(Q_ij) == 0)
  Qprime_ij[,which_pp] = 0
  
  # Solve for trophic level
  pseudoinvIminusQ = pseudoinverse( diag(nrow(Qprime_ij)) - Qprime_ij )
  l_i = rep(1,nrow(Qprime_ij)) %*% pseudoinvIminusQ
  return( l_i )
}
