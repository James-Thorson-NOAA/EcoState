
#' @title Calculate tracers, e.g., trophic level
#'
#' @description Calculate how a tracer propagates through consumption.
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
#' @param inverse_method whether to use pseudoinverse or standard inverse
#' @param which_primary Which taxa have no consumption and are therefore primary producers
#' @param tracer_i an indicator matrix specifying the traver value
#'
#' @details
#' Trophic level \eqn{y_i} for each predator \eqn{i} is defined as:
#'
#' \deqn{ \mathbf{y = l Q^* + 1} }
#'
#' where \eqn{\mathbf{Q*}} is the proportion consumption for each predator (column)
#' of different prey (rows).  We identify primary producers as any taxa with no
#' consumption (a column of 0s), and assign them as the first trophic level.
#'
#' More generically, a tracer might be used to track movement of biomass through
#' consumption.  For example, if we have a tracer \eqn{x_i} that is 1 for the 
#' base of the pelagic food chain, and 0 otherwise, then we can calculate 
#' the proportion of pelagic vs. nonpelagic biomass for each taxon:
#'
#' \deqn{ \mathbf{y = l Q^* + x} }
#'
#' This then allows us to separate alternative components of the foodweb.
#'
#' @importFrom corpcor pseudoinverse
#'
#' @export
compute_tracer <-
function( Q_ij,
          inverse_method = c("Penrose_moore", "Standard"),
          which_primary = which(colSums(Q_ij)==0),
          tracer_i = rep(1,nrow(Q_ij)) ){
  # Start
  inverse_method = match.arg(inverse_method)
  assertDouble( tracer_i, lower=0, upper=1, len=nrow(Q_ij), any.missing=FALSE )
  
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
  x_i = t(tracer_i) %*% inverse_IminusQ
  return( x_i )
}
