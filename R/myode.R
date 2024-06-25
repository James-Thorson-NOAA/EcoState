
#' @title 
#' Generic ODE solver
#'
#' @description 
#' Interface for \code{RTMBode::ode}, itself an interface for \code{deSolve::ode}
#'
#' @param f function in the differential equation y' = f(x, y);
#'        defined as a function R \times R^m \rightarrow R^m, where m is the number of equations.
#' @param a starting time for the interval to integrate
#' @param b ending time for the interval to integrate.
#' @param y0 starting values at time \code{a}
#' @param steps the number of steps from a to b.
#' @param Pars named list of parameters passed to f
#'
#' @details
#' Todo
#'
#' @export
myode = function( f, a, b, y0, n, Pars, method = "lsode" ){
  
  # Designed for dBdt, which outputs vector (not list)
  f2 = function( Time, State, Pars ){
    out = f( Time=Time, State=State, Pars=Pars )
    return(list(out))
  }
  
  #out = RTMBode::ode( func = f2,
  out = deSolve::ode( func = f2,
                times = seq(a,b,length=n+1), # n+1 to match rk4
                y = y0,
                method = method,
                parms = Pars, 
                atol = 1e-8, 
                rtol = 1e-8  )
  return(list(x=out[,1], y=out[,-1]))
}
