#' @title 
#' Adams-Bashford-Moulton for system of equations
#'
#' @description 
#' Third-order Adams-Bashford-Moulton predictor-corrector method.
#'
#' @inheritParams rk4sys
#'
#' @details
#' Combined Adams-Bashford and Adams-Moulton (or: multi-step) method of third order 
#' with corrections according to the predictor-corrector approach.
#' Copied from pracma under GPL-3, with small modifications to work with RTMB
#'
#' @export
abm3pc_sys <-
function (f, a, b, y0, n, Pars) {
  if (!is.numeric(n) || length(n) != 1 || n < 5){ 
    stop("Argument 'n' must be an integer greater or equal to 5.")
  }
  n <- floor(n)
  m <- length(y0)
  h <- (b - a)/n
  k <- h/12
  x <- seq(a, b, by = h)
  z <- y <- matrix(nrow=n + 1, ncol=m)
  z[1,] <- f(a, y0, Pars)
  y[1,] <- y0
  k1 <- h * z[1,]
  k2 <- h * f(a + h/2, y0 + k1/2, Pars)
  k3 <- h * f(a + 0.75 * h, y0 + 0.75 * k2, Pars)
  y[2,] <- y0 + (2 * k1 + 3 * k2 + 4 * k3)/9
  z[2,] <- f(x[2], y[2,], Pars)
  k1 <- h * z[2,]
  k2 <- h * f(x[2] + h/2, y[2,] + k1/2, Pars)
  k3 <- h * f(x[2] + 0.75 * h, y[2,] + 0.75 * k2, Pars)
  y[3,] <- y[2,] + (2 * k1 + 3 * k2 + 4 * k3)/9
  z[3,] <- f(x[2], y[2,], Pars)
  zz <- yy <- matrix(nrow=n+1, ncol=m)
  #errorest <- numeric(n)
  for (i in 3:n) {
    yy[i+1,] <- y[i,] + k * (23 * z[i,] - 16 * z[i-1,] + 5 * z[i-2,])
    zz[i+1,] <- f(x[i+1], yy[i+1,], Pars)
    y[i+1,] <- y[i,] + k * (5 * zz[i+1,] + 8 * z[i,] - z[i-1,])
    z[i+1,] <- f(x[i+1], y[i+1,], Pars)
    #errorest[i+1] <- -0.1 * (y[i+1] - yy[i+1])
  }
  #errorest <- sqrt(abs(errorest))
  return(list(x = x, y = y))
}
