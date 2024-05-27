#' @title 
#' Non-stiff (and stiff) ODE solvers
#'
#' @description 
#' Runge-Kutta (2, 3)-method with variable step size, resp
#'
#' @param f function in the differential equation y' = f(x, y);
#'        defined as a function R \times R^m \rightarrow R^m, where m is the number of equations.
#' @param a starting time for the interval to integrate
#' @param b ending time for the interval to integrate.
#' @param y0 starting values at time \code{a}
#' @param n Not used
#' @param Pars named list of parameters passed to f
#' @param rtol relative tolerance.
#' @param atol absolute tolerance.
#'
#' @details
#' Copied from pracma under GPL-3, with small modifications to work with RTMB
#'
#' @export
ode23 <- 
function( f, 
          a, 
          b, 
          y0,
          n,
          Pars, 
          rtol = 0.001, 
          atol = 1e-06){

  if (is.vector(y0)) {
    y0 <- as.matrix(y0)
  }
  else if (is.matrix(y0)) {
    if (ncol(y0) != 1) stop("Argument 'y0' must be a vector or single column matrix.")
  }
  eps <- .Machine$double.eps
  realmin <- 1e-100
  tdir <- sign(b - a)
  threshold <- atol/rtol
  hmax <- abs(0.1 * (b - a))
  t <- a
  tout <- t
  y <- y0
  yout <- t(y)
  s1 <- f(t, y, Pars)
  r <- max(abs(s1/max(abs(y), threshold))) + realmin
  h <- tdir * 0.8 * rtol^(1/3)/r
  while (t != b) {
    hmin <- 16 * eps * abs(t)
    if (abs(h) > hmax) {
      h <- tdir * hmax
    }
    else if (abs(h) < hmin) {
      h <- tdir * hmin
    }
    if (1.1 * abs(h) >= abs(b - t)) h <- b - t
    s2 <- f(t + h/2, y + h/2 * s1, Pars)
    s3 <- f(t + 3 * h/4, y + 3 * h/4 * s2, Pars)
    tnew <- t + h
    ynew <- y + h * (2 * s1 + 3 * s2 + 4 * s3)/9
    s4 <- f(tnew, ynew, Pars)
    e <- h * (-5 * s1 + 6 * s2 + 8 * s3 - 9 * s4)/72
    err <- max(abs(e/max(max(abs(y), abs(ynew)), threshold))) + realmin
    if (err <= rtol) {
      t <- tnew
      y <- ynew
      tout <- c(tout, t)
      yout <- rbind(yout, t(y))
      s1 <- s4
    }
    h <- h * min(5, 0.8 * (rtol/err)^(1/3))
    if (abs(h) <= hmin) {
      warning("Step size too small.")
      t <- b
    }
  }
  return(list(t = c(tout), y = yout))
}
