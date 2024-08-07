
#' @title Calculate tracers, e.g., trophic level
#'
#' @description Calculate how a tracer propagates through consumption.
#'
#' @inheritParams ecostate
#'
#' @param Q_ij Consumption of each prey i by predator j in units biomass.
#' @param inverse_method whether to use pseudoinverse or standard inverse
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
#' @export
compute_tracer <-
function( Q_ij,
          inverse_method = c("Penrose_moore", "Standard"),
          type_i,
          tracer_i = rep(1,nrow(Q_ij)) ){
  # Start
  inverse_method = match.arg(inverse_method)
  assertDouble( tracer_i, lower=0, upper=1, len=nrow(Q_ij), any.missing=FALSE )
  
  # Indicators 
  which_primary = which( type_i=="auto" )
  which_detritus = which( type_i=="detritus" )
  
  # Rescale consumption
  colsums = colSums(Q_ij)
  # diag(vector) doesn't work when length(vector)=1, so using explicit construction
  #B = diag(1/colsums)
  B = matrix(0, nrow=length(colsums), ncol=length(colsums))
  B[cbind(seq_along(colsums),seq_along(colsums))] = 1/colsums
  #inverse_denom = rep(1,nrow(Q_ij)) %*% t(1/colsums)
  Qprime_ij = Q_ij %*% B

  # Identify primary producers
  Qprime_ij[,c(which_primary,which_detritus)] = 0
  
  # Solve for trophic level
  if( inverse_method == "Penrose_moore" ){
    #inverse_IminusQ = corpcor::pseudoinverse( diag(nrow(Qprime_ij)) - Qprime_ij )
    inverse_IminusQ = ginv( diag(nrow(Qprime_ij)) - Qprime_ij )
  }else{
    inverse_IminusQ = solve( diag(nrow(Qprime_ij)) - Qprime_ij )
  }
  x_i = t(tracer_i) %*% inverse_IminusQ
  return( x_i )
}

#' @title Penrose-Moore pseudoinverse
#' @description Extend \code{MASS:ginv} to work with RTMB
#' @param X Matrix used to compute pseudoinverse
#' @importFrom MASS ginv
ginv <- RTMB::ADjoint(function(X) {
    n <- sqrt(length(X))
    dim(X) <- c(n,n)
    MASS::ginv(X)
}, function(X,Y,dY) {
    n <- sqrt(length(X))
    dim(Y) <- dim(dY) <- c(n,n)
    -t(Y)%*%dY%*%t(Y)
}, name="ginv")
