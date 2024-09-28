
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

#' @title Elementwise product of sparse and dense matrices
#' @description Calculate elementwise product of sparse and dense matrices
#' @param Msparse sparse matrix
#' @param Mdense dense matrix
#' @export
elementwise_product = function(Msparse, Mdense){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  if(!("p" %in% names(attributes(Msparse)))) browser()
  if(length(Msparse@x)>0){
    Mout = Msparse
    j = rep( seq_len(nrow(Msparse)), times=diff(Msparse@p) )
    Mout@x = Msparse@x + Mdense[cbind(Msparse@i+1, j)]
  }else{
    Mout = Matrix::sparseMatrix(i=1, j=1, x=0)
  }
  return(Mout)
}

#' @export
adsparse_to_matrix = function(x){
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  y = matrix(0, nrow=nrow(x), ncol=ncol(x))
  j = rep( seq_len(nrow(x)), times=diff(x@p) )
  y[cbind(x@i+1, j)] = x@x
  y
}

#' @title Dirichlet-multinomial
#' @description Allows data-weighting as parameter
#' @examples
#' library(RTMB)
#' prob = rep(0.1,10)
#' x = rmultinom( n=1, prob=prob, size=20 )[,1]
#' f = function( ln_theta ) ddirmult(x, prob, ln_theta)
#' f( 0 )
#' F = MakeTape(f, 0)
#' F$jacfun()(0)
#'
#' @export
ddirmult <-
function( x,
          prob,
          ln_theta,
          log=TRUE ){

  # Pre-processing
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  Ntotal = sum(x)
  p_exp = prob / sum(prob)
  p_obs = x / Ntotal
  dirichlet_Parm = exp(ln_theta) * Ntotal
  logres = 0.0

  # https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/main/SS_objfunc.tpl#L306-L314
  # https://github.com/James-Thorson/CCSRA/blob/master/inst/executables/CCSRA_v8.cpp#L237-L242
  # https://www.sciencedirect.com/science/article/pii/S0165783620303696

  # 1st term -- integration constant that could be dropped
  logres = logres + lgamma( Ntotal+1 )
  for( index in seq_along(x) ){
    logres = logres - lgamma( Ntotal*p_obs[index] + 1.0 )
  }

  # 2nd term in formula
  logres = logres + lgamma( dirichlet_Parm ) - lgamma( Ntotal+dirichlet_Parm )

  # Summation in 3rd term
  for( index in seq_along(x) ){
    logres = logres + lgamma( Ntotal*p_obs[index] + dirichlet_Parm*p_exp[index] )
    logres = logres - lgamma( dirichlet_Parm * p_exp[index] )
  }
  if(log){ return(logres) }else{ return(exp(logres)) }
}

combine_groups <-
function( taxa_list,
          B,
          PB,
          QB,
          DC ){

  # Aggregate rates
  PB_i = sapply( taxa_list, FUN=\(x){
               if(length(x)==1){
                 PB[match(x,names(PB))]
               }else{
                 weighted.mean(x=PB[match(x,names(PB))], w=B[match(x,names(B))], na.rm=TRUE)
               }
             } )
  QB_i = sapply( taxa_list, FUN=\(x){
               if(length(x)==1){
                 QB[match(x,names(PB))]
               }else{
                 weighted.mean(x=QB[match(x,names(QB))], w=B[match(x,names(B))])
               }
             } )
  B_i = sapply( taxa_list, FUN=\(x){
               sum(B[match(x,names(B))],na.rm=TRUE)
             } )
  names(PB_i) = names(QB_i) = names(taxa_list)

  # Aggregate diet matrix
  DC2 = t(sapply( taxa_list, FUN=\(x){
               colSums(DC[match(x,rownames(DC)),,drop=FALSE],na.rm=TRUE)
             } ))
  DC_ij = sapply( taxa_list, FUN=\(x){
               w = B[match(x,names(B))]
               rowSums(DC2[,match(x,colnames(DC2)),drop=FALSE] * outer(rep(1,nrow(DC2)),w) / sum(w), na.rm=TRUE )
             } )

  out = list(
    PB_i = PB_i,
    QB_i = QB_i,
    B_i = B_i,
    DC_ij = DC_ij
  )
  return(out)
}
