# cost_function ----------------------------------------------------------------
#' Cost function for optimization.
#'
#' \code{cost_function}
#'
#' @param y Factors to optimize.
#' @param rsm_output Second order model fit.
#' @param nfact Number of factor values.
#' @param xnam Formula names for \code{\link{rsm}}.
#'
#' @keywords internal
#' @export
cost_function = function(y, rsm_output, nfact, xnam) {
  tmp = vector(mode = "character", length = nfact)
  for (i in 1:nfact) {
    tmp[[i]] = paste0(xnam[i],"=y[",i,"]")
  }
  factors_coded = eval(parse(text = paste("data.frame(",paste(tmp, collapse = ","), ")")))
  abs(stats::predict(rsm_output, factors_coded))
}

# desirability -----------------------------------------------------------------
#' Quadratic desirability function.
#'
#' \code{desirability} is the cost function used for optimization.
#'
#' @param factors Factors (tuning variables).
#' @param Target Target value for the response.
#' @param a Parameter of the quadratic desirabilty function, defaults to a = 0.
#' @param b Parameter of the quadratic desirabilty function, defaults to b = 0.2.
#' @param c Parameter of the quadratic desirabilty function, defaults to c = 0.8.
#' @param Pred_min Minimum value which can be reached by the model within the factor ranges.
#' @param Pred_max Maximum value which can be reached by the model within the factor ranges.
#' @param nfact Number of factor values.
#' @param xnam Formula names for \code{\link{rsm}}.
#' @param tuner_rsm Second order model fit.
#'
#' @keywords internal
#' @export
desirability = function(factors, Target, w, a, b, c, nfact, xnam, tuner_rsm,
                        Pred_min, Pred_max) {

  y = cost_function(factors, tuner_rsm, nfact, xnam)

  if (y < Target) {
    P = Pred_min
  } else {
    P = Pred_max
  }

  Target = min(Pred_max, max(Pred_min, Target))

  X = (Target - y)/(Target - P)

  desirability = a + b*X + c*X^2 - 1

  if (is.na(desirability)) { desirability = 0 }  # e.g. when Target = P

  return(desirability)
}

# desirability_overall ---------------------------------------------------------
#' Overall desirability function.
#'
#' \code{desirability_overall} is the cost function for optimization with more
#' than one response variable. It is the weighted mean of the individual
#' desirability functions.
#'
#' @param factors Factors (tuning variables).
#' @param Target Vector of the target values for each response.
#' @param w Vector of weights for each the response.
#' @param nfact Number of factor values.
#' @param xnam Formula names for \code{\link{rsm}}.
#' @param tuner_rsm Second order model fit.
#' @param Pred_min Vector of minimum values which can be reached by the model
#' within the factor ranges for each response variable.
#' @param Pred_max Vector of maximum value which can be reached by the model
#' within the factor ranges for each response variable.
#'
#' @keywords internal
#' @export
desirability_overall = function(factors, Target, w, nfact, xnam, tuner_rsm,
                                Pred_min, Pred_max) {

  desir = numeric(length(w))
  for (i in 1:length(w)) {
    desir[i] = desirability(factors = factors, Target = Target[i],
                            w = w[i], a = 0, b = 0.2, c = 0.8, nfact = nfact,
                            xnam = xnam, tuner_rsm = tuner_rsm[[i]],
                            Pred_min = Pred_min[i], Pred_max = Pred_max[i])
  }

  desirability_overall = sum(w*desir)/sum(w)

  return(desirability_overall)
}

# GLPM axial potential --------------------------------------------------------------
#' Calculates the axial potential of gridless planar mirrors.
#'
#' \code{glpm_potential} calculates the axial potential of gridless planar mirrors.
#'
#' @param x Axial distance from the end electrode. x is normalized to the height
#' H of the mirror electrodes (same height for all electrodes).
#' @param L Vector of electrode lengths normalized with H ("length"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0).
#'
#' @return Axial potential at \code{x}.
#'
#' @references Yavor, M.I. et al. (2018), High performance gridless ion mirrors 
#' for multi-reflection time-of-flight and electrostatic trap mass analyzers, 
#' \emph{International Journal of Mass Spectrometry}, \strong{426},
#' 1-11, doi:10.1016/j.ijms.2018.01.009.
#' 
#' @keywords internal
#' @export
glpm_potential = function(x, L, V) {
  a = L
  b = L
  n = length(L)
  a[1] = 0
  for (i in 1:(n-1)) {
    a[i+1] = sum(L[1:i])
  }
  for (i in 1:n) {
    b[i] = sum(L[1:i])
  }
  
  Vx = 4*V[1]/pi*atan(exp(-pi*x)) 
  for (i in 1:n) {
    Vx = Vx + 2*V[i]/pi*(atan(exp(pi*(x-a[i]))) + atan(exp(pi*(x+a[i])))) -
      2*V[i]/pi*(atan(exp(pi*(x-b[i]))) + atan(exp(pi*(x+b[i]))))
  }
  return(Vx)
}

# GLPM inverse of axial potential ---------------------------------------------------
#' Calculates the inverse of the axial potential function.
#'
#' \code{glpm_potential_inv} calculates the position x given the axial potential (inverse
#' of axial potential function).
#'
#' The inverse of the axial potential function is used to calculate the turning 
#' point inside the mirror for a given kinetic energy.
#' 
#' @param y Axial potential.
#' @param L Vector of electrode lengths normalized with the height H of the 
#' mirror electrodes ("length"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#'
#' @return Position x (in units of H) where the axial potential = y.
#' 
#' @keywords internal
#' @export
glpm_potential_inv = function(y, L, V) {
  stats::uniroot((function(x,L,V) glpm_potential(x,L,V)-y), interval = c(0,5), 
                 L, V, tol = 1e-15)$root
}

# GLPM integrand for tof period calculation -----------------------------------------
#' Integrand for tof period calculation.
#'
#' \code{glpm_integrand} is the integrand for the time-of-flight period calculation.
#' 
#' @param x Axial distance from the end electrode, normalized to the height H of
#' the mirror electrodes.
#' @param E Potential energy at the turning point (normalized with the mean 
#' energy (K/K_0).
#' @param L Vector of electrode lengths normalized with H ("length"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#'
#' @return Axial potential at \code{x}.
#' 
#' @keywords internal
#' @export
glpm_integrand = function(x, E, L, V) {
  
  a = L
  b = L
  n = length(L)
  a[1] = 0
  for (i in 1:(n-1)) {
    a[i+1] = sum(L[1:i])
  }
  for (i in 1:n) {
    b[i] = sum(L[1:i])
  }
  Vx = 4*V[1]/pi*atan(exp(-pi*x)) 
  for (i in 1:n) {
    Vx = Vx + 2*V[i]/pi*(atan(exp(pi*(x-a[i]))) + atan(exp(pi*(x+a[i])))) -
      2*V[i]/pi*(atan(exp(pi*(x-b[i]))) + atan(exp(pi*(x+b[i]))))
  }
  Vx = 1/sqrt(E - Vx)
  
  return(Vx)
}

# GLPM tof period -------------------------------------------------------------------
#' Time-of-flight period calculation.
#'
#' \code{glpm_tofperiod} calculates the time-of-flight period inside the mirror.
#' 
#' This integrates the particle motion from the turning point x0 to the distance
#' x1. Because here we are only interested in relative time-of-flight deviations
#' all constant factors (e.g. \code{sqrt(m*amu/e/2)}) are omitted and the 
#' time-of-flight is in arbitrary units.
#' 
#' @param E Vector of potential energies at the turning point (normalized with
#' the mean energy (K/K_0).
#' @param x1 End distance of the particle, typically the time-of-flight focal 
#' point of the mirror (in units of H).
#' @param L Vector of electrode lengths normalized with the height H of the 
#' mirror electrodes ("length"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
glpm_tofperiod = function(E, x1, L, V) {
  x0 = sapply(E, glpm_potential_inv, L, V)
  tof = mapply(function(x0, x1, E) stats::integrate(glpm_integrand, x0, x1, E, L, V,
                                                    rel.tol = 1e-7)$value, x0, x1, E)
  return(tof)
}

# GLPM tof focal point --------------------------------------------------------------
#' Finds the time-of-flight focal point of the mirror.
#'
#' \code{glpm_find_x1} finds the time-of-flight focal point of the mirror.
#' 
#' @param L Vector of electrode lengths normalized with the height H of the 
#' mirror electrodes ("length"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#'
#' @return focal distance, measured from the back plane, normalized with H.
#' 
#' @keywords internal
#' @export
glpm_find_x1 = function(L, V) {
  E = seq(0.999, 1.001, length.out = 2)
  x1 = 3
  dx = 1
  repeat {
    tmp = glpm_tofperiod(E = E, x1 = x1, L,V)
    if (diff(tmp) < 0) {
      x1 = x1 - dx
      dx = dx/10
    }
    x1 = x1 + dx
    if (dx < 1e-6) break
  }
  return(x1)
}

# ZEIM variables in package environment ----------------------------------------
# https://stackoverflow.com/questions/12598242/global-variables-in-packages-in-r
pkg_env = new.env()
pkg_env$m = 31
pkg_env$X0 = gsl::bessel_zero_J0(1:pkg_env$m)
pkg_env$J1 = besselJ(pkg_env$X0,1)

# ZEIM axial potential --------------------------------------------------------------
#' Calculates the axial potential of cylindrical Zhang-Enke ion mirrors.
#'
#' \code{zeim_potential} calculates the axial potential of cylindrical 
#' three-element Zhang-Enke ion mirrors.
#'
#' @param x Axial distance from the entrance grid.
#' @param z1 Distance where electrode 1 starts.
#' @param z2 Distance where electrode 2 starts.
#' @param L Total length of the mirror, i.e. disance where end cap is.
#' @param V1 Voltage of electrode 1.
#' @param V2 Voltage of electrode 2.
#' @param R Radius of the cylinder.
#'
#' @return Axial potential at \code{x}.
#'
#' @references Yefchak, G.E. and Flory, C.A. (2002), Improved method for 
#' designing a cylindrical Zhang-Enke ion mirror, \emph{International Journal 
#' of Mass Spectrometry}, \strong{214}, 89-94, doi:10.1016/S1387-3806(01)00564-4.
#' 
#' @keywords internal
#' @export
zeim_potential = function(x, z1, z2, L, V1, V2, R) {
  
  # Note: x needs to be of length 1.
  I0 = besselI((1:pkg_env$m)*pi*R/L, 0)
  tmp = sinh(pkg_env$X0/R*x)/(pkg_env$X0*pkg_env$J1*sinh(pkg_env$X0/R*L))
  # Note: sum(tmp) does not converge well for x=L -> add half of last summand
  # i.e. (sum(tmp[1:m] + sum(tmp[2:m+1]))/2
  Vx = 2*V2*(sum(tmp[1:(pkg_env$m-1)], na.rm = TRUE) + tmp[pkg_env$m]/2) +
    2*V1/pi*sum((cos((1:pkg_env$m)*pi*z1/L)-cos((1:pkg_env$m)*pi*z2/L))/
                  ((1:pkg_env$m)*I0)*sin((1:pkg_env$m)*pi*x/L), na.rm = TRUE) +
    2*V2/pi*sum((cos((1:pkg_env$m)*pi*z2/L)-(-1)^(1:pkg_env$m))/
                  ((1:pkg_env$m)*I0)*sin((1:pkg_env$m)*pi*x/L), na.rm = TRUE)
  return(Vx)
}

# ZEIM inverse of axial potential ---------------------------------------------------
#' Calculates the inverse of the axial potential function.
#'
#' \code{zeim_potential_inv} calculates the position x given the axial potential (inverse
#' of axial potential function).
#'
#' The inverse of the axial potential function is used to calculate the turning 
#' point inside the mirror for a given kinetic energy.
#' 
#' @param x Axial distance from the entrance grid.
#' @param z1 Distance where electrode 1 starts.
#' @param z2 Distance where electrode 2 starts.
#' @param L Total length of the mirror, i.e. disance where end cap is.
#' @param V1 Voltage of electrode 1.
#' @param V2 Voltage of electrode 2.
#' @param R Radius of the cylinder.
#'
#' @return Position x where the axial potential = y.
#' 
#' @keywords internal
#' @export
zeim_potential_inv = function(y, z1, z2, L, V1, V2, R) {
  stats::uniroot((function(x,z1,z2,L,V1,V2,R) zeim_potential(x,z1,z2,L,V1,V2,R)-y), 
                 interval = c(0,L), z1, z2, L, V1, V2, R, tol = 1e-15)$root
}

# ZEIM integrand for tof period calculation -----------------------------------------
#' Integrand for tof period calculation.
#'
#' \code{zeim_integrand} is the integrand for the time-of-flight period calculation.
#' 
#' @param x Axial distance from the entrance grid.
#' @param E Potential energy at the turning point.
#' @param z1 Distance where electrode 1 starts.
#' @param z2 Distance where electrode 2 starts.
#' @param L Total length of the mirror, i.e. disance where end cap is.
#' @param V1 Voltage of electrode 1.
#' @param V2 Voltage of electrode 2.
#' @param R Radius of the cylinder.
#'
#' @return Axial potential at \code{x}.
#' 
#' @keywords internal
#' @export
zeim_integrand = function(x, E, z1, z2, L, V1, V2, R) {
  
  I0 = besselI((1:pkg_env$m)*pi*R/L, 0)
  a1 = 0
  a2 = 0
  a3 = 0
  for (i in 1:(pkg_env$m-1)) {
    a1 = a1 + sinh(pkg_env$X0[i]/R*x)/(pkg_env$X0[i]*pkg_env$J1[i]*sinh(pkg_env$X0[i]/R*L))
  }
  a1 = a1 + sinh(pkg_env$X0[pkg_env$m]/R*x)/
    (pkg_env$X0[pkg_env$m]*pkg_env$J1[pkg_env$m]*sinh(pkg_env$X0[pkg_env$m]/R*L))/2
  for (i in 1:pkg_env$m) {
    a2 = a2 + (cos(i*pi*z1/L)-cos(i*pi*z2/L))/(i*I0[i])*sin(i*pi*x/L)
  }
  for (i in 1:pkg_env$m) {
    a3 = a3 + (cos(i*pi*z2/L)-(-1)^i)/(i*I0[i])*sin(i*pi*x/L)
  }
  Vx = 2*V2*a1 + 2*V1/pi*a2 + 2*V2/pi*a3
  
  Vx = 1/sqrt(E - Vx)
  
  return(Vx)
}

# ZEIM tof period -------------------------------------------------------------------
#' Time-of-flight period calculation.
#'
#' \code{zeim_tofperiod} calculates the time-of-flight period inside the mirror.
#' 
#' This integrates the particle motion from the entrance grid to the turning 
#' point x0. Because here we are only interested in relative time-of-flight 
#' deviations all constant factors (e.g. \code{sqrt(m*amu/e/2)}) are omitted 
#' and the time-of-flight is in arbitrary units.
#' 
#' @param E Potential energy at the turning point.
#' @param z1 Distance where electrode 1 starts.
#' @param z2 Distance where electrode 2 starts.
#' @param L Total length of the mirror, i.e. disance where end cap is.
#' @param V1 Voltage of electrode 1.
#' @param V2 Voltage of electrode 2.
#' @param R Radius of the cylinder.
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
zeim_tofperiod = function(E, z1, z2, L, V1, V2, R) {
  x0 = sapply(E, zeim_potential_inv, z1, z2, L, V1, V2, R)
  tof = mapply(function(x0, E) stats::integrate(zeim_integrand, 0, x0, E, z1, 
                                                z2, L, V1, V2, R, rel.tol = 1e-6)$value, x0, E)
  return(tof)
}

# ZEIM tof focus point --------------------------------------------------------------
#' Finds the time-of-flight focal point of the mirror.
#'
#' \code{zeim_find_x1} finds the time-of-flight focal point of the mirror.
#' 
#' @param z1 Distance where electrode 1 starts.
#' @param z2 Distance where electrode 2 starts.
#' @param L Total length of the mirror, i.e. disance where end cap is.
#' @param V1 Voltage of electrode 1.
#' @param V2 Voltage of electrode 2.
#' @param R Radius of the cylinder.
#'
#' @return focal point (measured from the entrance grid).
#' 
#' @keywords internal
#' @export
zeim_find_x1 = function(z1, z2, L, V1, V2, R) {
  E = seq(0.99, 1.01, length.out = 2)
  x1 = L
  dx = 10
  repeat {
    tmp = zeim_tofperiod(E = E, z1, z2, L, V1, V2, R) + x1*1/sqrt(E)
    if (diff(tmp) < 0) {
      x1 = x1 - dx
      dx = dx/10
    }
    x1 = x1 + dx
    if (dx < 1e-5) break
  }
  return(x1)
}
# PIM basic potential function -------------------------------------------------
U0 = function(V, z) {
  2*V/pi*atan(exp(pi*z))
}

# PIM axial potential --------------------------------------------------------------
#' Calculates the axial potential of a planar multi-element ion mirror.
#'
#' \code{pim_potential} calculates the axial potential of a planar 
#' multi-element ion mirror with entrance grid.
#'
#' @param x Axial distance from the entrance grid (normalized with H).
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} is
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of normalized electrode potentials. The entrance grid and the 
#' first electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized to the height of the mirror electrodes H 
#' (which is the same for all electrodes) and the potentials are normalized to
#' the mean kinetic energy ("voltage"/K_0) of the (singly-charged) ions.
#'
#' @return Axial potential at \code{x}.
#'
#' @references Yavor, M.I. et al. (2018), High performance gridless ion mirrors 
#' for multi-reflection time-of-flight and electrostatic trap mass analyzers, 
#' \emph{International Journal of Mass Spectrometry}, \strong{426},
#' 1-11, doi:10.1016/j.ijms.2018.01.009.
#' 
#' @keywords internal
#' @export
pim_potential = function(x, z, V) {
  Vx = rep(0, length(x))
  n = length(z)
  if (n<2) stop("z and V need to be vectors with at least two elements.")
  for (i in 1:(n-1)) {
    # Vx = Vx + U0(V[i+1], x-z[i]) - U0(V[i+1], x-z[i+1]) -
    # U0(V[i+1], -x-z[i]) + U0(V[i+1], -x-z[i+1]) -
    #   U0(V[i+1], -x+2*z[n]-z[i]) + U0(V[i+1], -x+2*z[n]-z[i+1]) +
    #   U0(V[i+1], x+2*z[n]-z[i]) - U0(V[i+1], x+2*z[n]-z[i+1]) +
    #   U0(V[i+1], x-2*z[n]-z[i]) - U0(V[i+1], x-2*z[n]-z[i+1]) -
    #   U0(V[i+1], -x-2*z[n]-z[i]) + U0(V[i+1], -x-2*z[n]-z[i+1]) -
    #   U0(V[i+1], -x+4*z[n]-z[i]) + U0(V[i+1], -x+4*z[n]-z[i+1]) +
    #   U0(V[i+1], x+4*z[n]-z[i]) - U0(V[i+1], x+4*z[n]-z[i+1]) +
    #   U0(V[i+1], x-4*z[n]-z[i]) - U0(V[i+1], x-4*z[n]-z[i+1]) -
    #   U0(V[i+1], -x-4*z[n]-z[i]) + U0(V[i+1], -x-4*z[n]-z[i+1]) -
    #   U0(V[i+1], -x+6*z[n]-z[i]) + U0(V[i+1], -x+6*z[n]-z[i+1])
    Vx = Vx + U0(V[i+1], x-z[i]) - U0(V[i+1], x-z[i+1])
    for (j in 1:4) {
      Vx = Vx - (U0(V[i+1],-x+2*j*z[n]-z[i]) - U0(V[i+1],-x+2*j*z[n]-z[i+1]) +
                   U0(V[i+1],-x-2*(j-1)*z[n]-z[i]) - U0(V[i+1],-x-2*(j-1)*z[n]-z[i+1])) +
        (U0(V[i+1],x+2*j*z[n]-z[i]) - U0(V[i+1],x+2*j*z[n]-z[i+1]) +
           U0(V[i+1],x-2*j*z[n]-z[i]) - U0(V[i+1],x-2*j*z[n]-z[i+1]))
    }
  }
  Vx = Vx + 2*(U0(V[n], x-z[n]) - U0(V[n], -x-z[n]))  # end cap
  
  return(Vx)
}

# PIM inverse of axial potential ---------------------------------------------------
#' Calculates the inverse of the axial potential function.
#'
#' \code{pim_potential_inv} calculates the position x given the axial potential (inverse
#' of axial potential function).
#'
#' The inverse of the axial potential function is used to calculate the turning 
#' point inside the mirror for a given kinetic energy.
#' 
#' @param x Axial distance from the entrance grid (normalized with H).
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of normalized electrode potentials. The entrance grid and the 
#' first electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized to the height of the mirror electrodes H 
#' (which is the same for all electrodes) and the potentials are normalized to
#' the mean kinetic energy ("voltage"/K_0) of the (singly-charged) ions.
#'
#' @return Position x where the axial potential = y.
#' 
#' @keywords internal
#' @export
pim_potential_inv = function(y, z, V) {
  stats::uniroot((function(x,z,V) pim_potential(x,z,V)-y), 
                 interval = c(0,z[length(z)]), z, V, tol = 1e-15)$root
}

# PIM integrand for tof period calculation -----------------------------------------
#' Integrand for tof period calculation.
#'
#' \code{pim_integrand} is the integrand for the time-of-flight period calculation.
#' 
#' @param x Axial distance from the entrance grid (normalized with H).
#' @param E Potential energy at the turning point (normalized with H).
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of normalized electrode potentials. The entrance grid and the 
#' first electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized to the height of the mirror electrodes H 
#' (which is the same for all electrodes) and the potentials are normalized to
#' the mean kinetic energy ("voltage"/K_0) of the (singly-charged) ions.
#'
#' @return Axial potential at \code{x}.
#' 
#' @keywords internal
#' @export
pim_integrand = function(x, E, z, V) {

  Vx = 1/sqrt(E - pim_potential(x, z, V))
  
  return(Vx)
}

# PIM tof period -------------------------------------------------------------------
#' Time-of-flight period calculation.
#'
#' \code{pim_tofperiod} calculates the time-of-flight period inside the mirror.
#' 
#' This integrates the particle motion from the entrance grid to the turning 
#' point x0. Because here we are only interested in relative time-of-flight 
#' deviations all constant factors (e.g. \code{sqrt(m*amu/e/2)}) are omitted 
#' and the time-of-flight is in arbitrary units.
#' 
#' @param E Potential energy at the turning point.
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of normalized electrode potentials. The entrance grid and the 
#' first electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized to the height of the mirror electrodes H 
#' (which is the same for all electrodes) and the potentials are normalized to
#' the mean kinetic energy ("voltage"/K_0) of the (singly-charged) ions.
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
pim_tofperiod = function(E, z, V) {
  x0 = sapply(E, pim_potential_inv, z, V)
  tof = mapply(function(x0, E) stats::integrate(pim_integrand, 0, x0, E, z, V,
                                                rel.tol = 1e-6)$value, x0, E)
  return(tof)
}

# PIM total tof -------------------------------------------------------------------
#' Time-of-flight period calculation.
#'
#' \code{pim_totaltof} calculates the total time-of-flight including field-free
#' space and linear stage of the mirror.
#' 
#' Because here we are only interested in relative time-of-flight 
#' deviations all constant factors (e.g. \code{sqrt(m*amu/e/2)}) are omitted 
#' and the time-of-flight is in arbitrary units.
#' 
#' @param E Potential energy at the turning point.
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of normalized electrode potentials. The entrance grid and the 
#' first electrode are assumed to be grounded (\code{V[1]=0}).
#' @param x1 Focal distance of mirror (normalized with H), measured from the 
#' entrance grid).
#' @param d5 Distance of first linear stage of mirror (normalized with H).
#' @param u5 Normalized potential in first linear stage of mirror.
#' 
#' All distances are normalized to the height of the mirror electrodes H 
#' (which is the same for all electrodes). The potentials are normalized to
#' the mean kinetic energy ("voltage"/K_0) of the (singly-charged) ions with
#' \code{u5=0}. With \code{u5!=0}, the normalized kinetic energy is \code{1+u5}.
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
pim_totaltof = function(E, z, V, x1, d5 = 0, u5 = 0) {
  if (u5 != 0) {
    tof = pim_tofperiod(E, z, V) + x1*1/sqrt(E+u5) + 2*d5/u5*(sqrt(E+u5)-sqrt(E))
  } else {
    tof = pim_tofperiod(E, z, V) + x1*1/sqrt(E)
  }
  # note: this is actually only half of the tof: from x1 to the return point in
  # the mirror.
  return(tof)
}

# PIM tof focal point --------------------------------------------------------------
#' Finds the time-of-flight focal point of the mirror.
#'
#' \code{pim_find_x1} finds the time-of-flight focal point of the mirror.
#' 
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of normalized electrode potentials. The entrance grid and the 
#' first electrode are assumed to be grounded (\code{V[1]=0}).
#' @param d5 Distance of first linear stage of mirror (normalized with H).
#' @param u5 Normalized potential in first linear stage of mirror.
#' 
#' All distances are normalized to the height of the mirror electrodes H 
#' (which is the same for all electrodes). The potentials are normalized to
#' the mean kinetic energy ("voltage"/K_0) of the (singly-charged) ions with
#' \code{u5=0}. With \code{u5!=0}, the normalized kinetic energy is \code{1+u5}.
#'
#' @return focal distance, measured from the entrance grid.
#' 
#' @keywords internal
#' @export
pim_find_x1 = function(z, V, d5 = 0, u5 = 0) {
  E = seq(0.99, 1.01, length.out = 2)
  x1 = 0
  dx = 10
  repeat {
    tmp = pim_totaltof(E, z, V, x1, d5, u5)
    if (diff(tmp) < 0) {
      x1 = x1 - dx
      dx = dx/10
    }
    x1 = x1 + dx
    if (dx < 1e-5) break
  }
  if (x1<0) stop("x1 < 0")
  return(x1)
}

# SNLM correction integrand ----------------------------------------------------
snlm_int = function(y, a, b, u1, u3, u5, u6, d1, d3, d4, d5, d6, x0) {
  k = x0/d1
  f1 = 2*d1/sqrt(u1)*(sqrt(k)) +
       2*d3/u3*(sqrt(k*u1+u3)-sqrt(k*u1)) +
       d4/sqrt(k*u1+u3) +
       4*d5/u5*(sqrt(k*u1+u3)-sqrt(k*u1+u3-u5)) +
       4*d6/u6*sqrt(k*u1+u3-u5)
  x = (a-y)*d1/u1
  k = (x0-x)/d1
  f2 = 2*d1/sqrt(u1)*(sqrt(k)) +
       2*d3/u3*(sqrt(k*u1+u3)-sqrt(k*u1)) +
       d4/sqrt(k*u1+u3) +
       4*d5/u5*(sqrt(k*u1+u3)-sqrt(k*u1+u3-u5)) +
       4*d6/u6*sqrt(k*u1+u3-u5)
  f = (f1 - f2)/sqrt(b-y)
  return(f)
}

# SNLM potential ---------------------------------------------------------------
#' Calculates the potential of a Shimadzu non-linear mirror.
#'
#' \code{snlm_potential} calculates the potential of a Shimadzu non-linear 
#' mirror, which is a linear two-stage ion mirror with a non-linear correction
#' potential in the second stage.
#'
#' @param d1 Distance of first stage of extractor, from push to pull (m).
#' @param d3 Distance of second stage of extractor, from pull to drift (m).
#' @param d4 Distance of of drift space, from extractor to reflector and from
#' reflector to detector (m).
#' @param d5 Distance of first stage of reflector (m).
#' @param d6 Distance of second stage of reflector (m).
#' @param shift Reference energy shift inside extractor (m).
#' @param drift Drift tube voltage (V).
#' @param pulse Extraction pulse voltage (V).
#'
#' @return A data.frame with the corrected total distance, the corresponding 
#' potential and the distance correction in the second stage of the mirror.
#'
#' @references U.S. Patent No. 8,772,708 B2 (filed Dec. 20, 2011).
#' 
#' @keywords internal
#' @export
snlm_potential = function(d1, d3, d4, d5, d6, shift, drift, pulse) {
  
  u1 = 2*pulse
  u3 = -pulse - drift
  x0 = d1 - 0.001 - shift
  k0 = x0/d1
  p0 = k0 + u3/u1
  a = d1/u1*k0^(-1/2) + d3/u3*(p0^(-1/2)-k0^(-1/2)) - d4/u1/2*p0^(-3/2)
  b = d1/u1/2*k0^(-3/2) + d3/u3/2*(p0^(-3/2)-k0^(-3/2)) - d4/u1*3/4*p0^(-5/2)
  u5 = (a - 2*p0*b + 2*d5/u1*p0^(-3/2))/(-2*b/u1)
  u6 = (-2*d6*(p0-u5/u1)^(-1/2))/(a + 2*d5/u5*(p0^(-1/2)-(p0-u5/u1)^(-1/2)))
  
  E0 = u1*x0/d1+u3
  
  # double-exponential transformation for fast numerical itegration
  a = E0  # lower integration boundary
  h = 0.1  # step size   
  l = seq(0+h/2, 2.5-h/2, h)
  phi = tanh(pi/2*sinh(l))
  dphi = pi/2*cosh(l)/(cosh(pi/2*sinh(l))^2)
  
  U = seq(E0+1, max(u1+u3, u5+u6), 1)
  if ((u5+u6)<(u1+u3)) { # reflector too short
    d6plus = d6*((u1+u3) - (u5+u6))/u6
    print(paste("reflector needs to be longer by at least", d6plus, "m"))
  }
  xc = rep(0, length(U))
  
  for (i in 1:length(U)) {
    b = U[i]  # upper integration boundary
    c = (b-a)/2  # variable transformation
    d = (a+b)/2  # variable transformation
    xc[i] = 1/2/pi*c*h*
      sum((snlm_int(c*phi+d,a,b,u1,u3,u5,u6,d1,d3,d4,d5,d6,x0) +
             snlm_int(-c*phi+d,a,b,u1,u3,u5,u6,d1,d3,d4,d5,d6,x0))*dphi)
  }
  
  xa = d6*(U-u5)/u6
  x = xa + xc
  U0 = u6/d6*xa
  ret = data.frame(x = x, U = U0, xc = xc)
  return(ret)
}