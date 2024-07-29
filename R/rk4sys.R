#' @title 
#' Classical Runge-Kutta for system of equations
#'
#' @description 
#' Classical Runge-Kutta of order 4.
#'
#' @param f function in the differential equation \eqn{y' = f(x, y)};
#'        defined as a function \eqn{R \times R^m \rightarrow R^m}, where \eqn{m} is the number of equations.
#' @param a starting time for the interval to integrate
#' @param b ending time for the interval to integrate.
#' @param y0 starting values at time \code{a}
#' @param steps the number of steps from a to b.
#' @param Pars named list of parameters passed to f
#'
#' @details
#' Classical Runge-Kutta of order 4 for (systems of) ordinary differential equations with fixed step size.
#' Copied from pracma under GPL-3, with small modifications to work with RTMB
#'
#' @export
rk4sys <-
function (f, a, b, y0, n, Pars) {
  m <- length(y0)
  h <- (b - a)/n
  x <- seq(a + h, b, by = h)
  y <- matrix(0, nrow = n, ncol = m)
  k1 <- h * f(a, y0, Pars)
  k2 <- h * f(a + h/2, y0 + k1/2, Pars)
  k3 <- h * f(a + h/2, y0 + k2/2, Pars)
  k4 <- h * f(a + h, y0 + k3, Pars)
  y[1, ] <- y0 + k1/6 + k2/3 + k3/3 + k4/6
  for (i in seq_len(n-1)) {
    k1 <- h * f(x[i], y[i, ], Pars)
    k2 <- h * f(x[i] + h/2, y[i, ] + k1/2, Pars)
    k3 <- h * f(x[i] + h/2, y[i, ] + k2/2, Pars)
    k4 <- h * f(x[i] + h, y[i, ] + k3, Pars)
    y[i + 1, ] <- y[i, ] + k1/6 + k2/3 + k3/3 + k4/6
  }
  x <- c(a, x)
  y <- rbind(y0, y)
  return(list(x = x, y = y))
}
