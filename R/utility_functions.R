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

# plot_coeffs ------------------------------------------------------------------
#' Plots the fit coefficients of the response surface model.
#'
#' \code{plot_coeffs} generates a plot with the fit coefficients of the response
#' surface model with 95% confidence interval.
#'
#' @param tuner_rsm Second order model fit.
#' @param ylab y-axis label.
#' @param factors Factors.
#'
#' @keywords internal
#' @export
plot_coeffs = function(tuner_rsm, ylab = "", factors) {
  m = length(tuner_rsm$coeff)
  coeffs = tuner_rsm$coeff[2:m]
  stderror = summary(tuner_rsm)$coeff[2:m,2]
  confinterval = stats::confint(tuner_rsm)[2:m,]
  ylim = 1.2*range(coeffs)
  labels = sapply(strsplit(names(coeffs), ")", fixed = TRUE), "[", 2)
  nfact = length(factors$Name)
  for (i in 1:nfact) {
    labels = gsub(paste("x" ,i, sep = ""), factors$Name[i], labels)
  }
  mp = graphics::barplot(coeffs, axes=FALSE, axisnames=FALSE, ylim = ylim,
               main="Coefficients", xlab="", ylab=ylab,
               col = c(rep("blue", nfact), rep("darkgreen", sum(1:(nfact-1))), rep("red", nfact)))
  graphics::axis(1, labels = labels, at = mp, cex.axis = 0.6, las = 2)
  graphics::axis(2)
  graphics::box()
  # error bars
  graphics::segments(mp, confinterval[,1], mp, confinterval[,2])
  graphics::segments(mp - 0.2, confinterval[,1], mp + 0.2, confinterval[,1])
  graphics::segments(mp - 0.2, confinterval[,2], mp + 0.2, confinterval[,2])
}

# plot_results -----------------------------------------------------------------
#' Plots results.
#'
#' \code{plot_results} plots results.
#'
#' @param pltly Plotly object.
#' @param k Step number.
#' @param n_steps Total number of steps.
#' @param resultdir Results directory (used for legend text).
#' @param result Results (Resolving power and sensitivity).
#' @param runs Box-Behnken design data.
#' @param bestpoint_run Bestpoint control values.
#' @param bestpoint_result Verified bestpoint result.
#' @param bestpoint_predicted Model prediction of bestpoint result.
#' @param xylabel Axis labels.
#'
#' @keywords internal
#' @export
plot_results = function(pltly, k, n_steps, resultdir, result, runs, bestpoint_run,
                        bestpoint_result, bestpoint_predicted, xylabel) {

  colors = gplots::rich.colors(n_steps)

  legendtext = resultdir

  len = length(result[,1])
  
  # generate marker text strings
  tmp1 = list()
  for (i in 1:length(runs)) {
    tmp1[[i]] = paste0(names(runs)[i], ": ", signif(runs[,i], 12))
  }
  tmp2 = vector(mode = "character", length = length(runs[,1]))
  for (i in 1:length(runs[,1])) {
    tmp2[i] = paste0(lapply(tmp1, "[", i), collapse = "<br>")
  }
  tmp3 = list()
  for (i in 1:length(bestpoint_run)) {
    tmp3[[i]] = paste0(names(bestpoint_run)[i], ": ", signif(bestpoint_run[i], 12))
  }
  tmp4 = paste0(lapply(tmp3, "[", 1), collapse = "<br>")
  # add markers
  pltly = plotly::add_markers(pltly, data = result[1:(len-1),], x = ~res, y = ~sens,
                      marker = list(color = colors[k], symbol = "x"),
                      name = legendtext,
                      text = eval(tmp2[1:(len-1)]))
  pltly = plotly::add_markers(pltly, data = result[len,], x = ~res, y = ~sens,
                      marker = list(color = colors[k], symbol = "triangle-up"),
                      name = legendtext, showlegend = FALSE,
                      text = eval(tmp2[len]))  # center point
  pltly = plotly::add_markers(pltly, data = bestpoint_result, x = ~res, y = ~sens,
                      marker = list(color = colors[k], symbol = "o", size = 12),
                      name = legendtext, showlegend = FALSE,
                      text = eval(tmp4))  # best point measured
  pltly = plotly::add_markers(pltly, data = bestpoint_predicted, x = ~res, y = ~sens,
                      marker = list(color = colors[k], symbol = "circle-open", size = 12),
                      name = legendtext, showlegend = FALSE,
                      text = eval(tmp4))  # best point predicted

  # plot
  pltly = plotly::layout(pltly, title = "", xaxis = list(title = xylabel[1]),
                 yaxis = list(title = xylabel[2]))
  print(pltly)
  return(pltly)
}

# GLPM axial potential --------------------------------------------------------------
#' Calculates the axial potential of gridless planar mirrors.
#'
#' \code{glpm_potential} calculates the axial potential of gridless planar mirrors.
#'
#' @param x Axial distance from the end electrode.
#' @param L Vector of electrode lengths normalized with H ("lenght"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#' @param H Height of the mirror electrodes (the same for all electrodes).
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
glpm_potential = function(x, L, V, H) {
  a = L
  b = L
  n = length(L)
  a[1] = 0
  for (i in 1:(n-1)) {
    a[i+1] = sum(L[1:i])*H
  }
  for (i in 1:n) {
    b[i] = sum(L[1:i])*H
  }
  
  Vx = 4*V[1]/pi*atan(exp(-pi*x/H)) 
  for (i in 1:n) {
    Vx = Vx + 2*V[i]/pi*(atan(exp(pi*(x-a[i])/H)) + atan(exp(pi*(x+a[i])/H))) -
      2*V[i]/pi*(atan(exp(pi*(x-b[i])/H)) + atan(exp(pi*(x+b[i])/H)))
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
#' @param L Vector of electrode lengths normalized with H ("lenght"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#' @param H Height of the mirror electrodes (the same for all electrodes).
#'
#' @return Position x where the axial potential = y.
#' 
#' @keywords internal
#' @export
glpm_potential_inv = function(y, L, V, H) {
  stats::uniroot((function(x,L,V,H) glpm_potential(x,L,V,H)-y), interval = c(0,100), 
                 L, V, H, tol = 1e-15)$root
}

# GLPM integrand for tof period calculation -----------------------------------------
#' Integrand for tof period calculation.
#'
#' \code{glpm_integrand} is the integrand for the time-of-flight period calculation.
#' 
#' @param x Axial distance from the end electrode.
#' @param E Potential energy at the turning point (normalized with the mean 
#' energy (K/K_0).
#' @param L Vector of electrode lengths normalized with H ("lenght"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#' @param H Height of the mirror electrodes (the same for all electrodes).
#'
#' @return Axial potential at \code{x}.
#' 
#' @keywords internal
#' @export
glpm_integrand = function(x, E, L, V, H) {
  
  a = L
  b = L
  n = length(L)
  a[1] = 0
  for (i in 1:(n-1)) {
    a[i+1] = sum(L[1:i])*H
  }
  for (i in 1:n) {
    b[i] = sum(L[1:i])*H
  }
  Vx = 4*V[1]/pi*atan(exp(-pi*x/H)) 
  for (i in 1:n) {
    Vx = Vx + 2*V[i]/pi*(atan(exp(pi*(x-a[i])/H)) + atan(exp(pi*(x+a[i])/H))) -
      2*V[i]/pi*(atan(exp(pi*(x-b[i])/H)) + atan(exp(pi*(x+b[i])/H)))
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
#' all constant factors (e.g. \code{2*sqrt(2*m*amu/e)}) are omitted and the 
#' time-of-flight is in arbitrary units.
#' 
#' @param E Vector of potential energies at the turning point (normalized with
#' the mean energy (K/K_0).
#' @param x1 End distance of the particle, typically the time-of-flight focal 
#' point of the mirror.
#' @param L Vector of electrode lengths normalized with H ("lenght"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#' @param H Height of the mirror electrodes (the same for all electrodes).
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
glpm_tofperiod = function(E, x1, L, V, H) {
  x0 = sapply(E, glpm_potential_inv, L, V, H)
  tof = mapply(function(x0, x1, E) stats::integrate(glpm_integrand, x0, x1, E, L, V, H,
                                                    rel.tol = 1e-6)$value, x0, x1, E)
  return(tof)
}

# GLPM tof focus point --------------------------------------------------------------
#' Finds the time-of-flight focal point of the mirror.
#'
#' \code{glpm_find_x1} finds the time-of-flight focal point of the mirror.
#' 
#' @param L Vector of electrode lengths normalized with H ("lenght"/\code{H}).
#' @param V Vector of voltages normalized with the mean energy ("voltage"/K_0)
#' @param H Height of the mirror electrodes (the same for all electrodes).
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
glpm_find_x1 = function(L, V, H) {
  E = seq(0.99, 1.01, length.out = 2)
  x1 = 3*H
  dx = 10
  repeat {
    tmp = glpm_tofperiod(E = E, x1 = x1, L,V,H)
    if (diff(tmp) < 0) {
      x1 = x1 - dx
      dx = dx/10
    }
    x1 = x1 + dx
    if (dx < 1e-5) break
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
#' Calculates the axial potential of cylindrical Zhang–Enke ion mirrors.
#'
#' \code{zeim_potential} calculates the axial potential of cylindrical 
#' three-element Zhang–Enke ion mirrors.
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
#' designing a cylindrical Zhang–Enke ion mirror, \emph{International Journal 
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
#' deviations all constant factors (e.g. \code{2*sqrt(2*m*amu/e)}) are omitted 
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
#' @return time-of-flight.
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
#' Calculates the axial potential of a planar three-element ion mirror.
#'
#' \code{pim_potential} calculates the axial potential of a planar 
#' three-element ion mirror.
#'
#' @param x Axial distance from the entrance grid (normalized with H).
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of electrode voltages. The entrance grid and the first
#' electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized with the height of the mirror electrodes H 
#' (which is the same for all electrodes).
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
#' @param V Vector of electrode voltages. The entrance grid and the first
#' electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized with the height of the mirror electrodes H 
#' (which is the same for all electrodes).
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
#' @param V Vector of electrode voltages. The entrance grid and the first
#' electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized with the height of the mirror electrodes H 
#' (which is the same for all electrodes).
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
#' deviations all constant factors (e.g. \code{2*sqrt(2*m*amu/e)}) are omitted 
#' and the time-of-flight is in arbitrary units.
#' 
#' @param E Potential energy at the turning point.
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of electrode voltages. The entrance grid and the first
#' electrode are assumed to be grounded (\code{V[1]=0}).
#' 
#' All distances are normalized with the height of the mirror electrodes H 
#' (which is the same for all electrodes).
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
#' deviations all constant factors (e.g. \code{2*sqrt(2*m*amu/e)}) are omitted 
#' and the time-of-flight is in arbitrary units.
#' 
#' @param E Potential energy at the turning point.
#' @param z Vector of electrode distances (normalized with H). \code{z[1]} it
#' the end of the entrance grid electrode and \code{z[length(z)]} is the total
#' length of the mirror.
#' @param V Vector of electrode voltages. The entrance grid and the first
#' electrode are assumed to be grounded (\code{V[1]=0}).
#' @param x1 Focal distance of mirror (normalized with H).
#' @param d5 Distance of first linear stage of mirror (normalized with H).
#' @param u5 Potential in first linear stage of mirror.
#' 
#' All distances are normalized with the height of the mirror electrodes H 
#' (which is the same for all electrodes).
#'
#' @return time-of-flight.
#' 
#' @keywords internal
#' @export
pim_totaltof = function(E, z, V, x1, d5 = 0, u5 = 0) {
  tof = pim_tofperiod(E, z, V) + x1*1/sqrt(E+u5) + 4*d5/u5*(sqrt(E+u5)-sqrt(E))
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
#' @param V Vector of electrode voltages. The entrance grid and the first
#' electrode are assumed to be grounded (\code{V[1]=0}).
#' @param d5 Distance of first linear stage of mirror (normalized with H).
#' @param u5 Potential in first linear stage of mirror.
#' 
#' All distances are normalized with the height of the mirror electrodes H 
#' (which is the same for all electrodes).
#'
#' @return focal point.
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