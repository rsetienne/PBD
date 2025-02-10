pbd_loglik_rhs_cpp = function(t,x,pars)
{
   b = pars[[1]](t,as.numeric(pars[5:(length(pars)-1)]))
   mu_g = pars[[2]](t,as.numeric(pars[5:(length(pars)-1)]))
   la2 = pars[[3]](t,as.numeric(pars[5:(length(pars)-1)]))
   mu_i = pars[[4]](t,as.numeric(pars[5:(length(pars)-1)]))
   nu2 = la2 + mu_i
   H = x[1]
   p_i = x[2]
   q_i = x[3]
   q_g = x[4]
   dH = -H * b * (1 - p_i)
   dp_i = -(nu2 + b) * p_i + la2 * q_g + mu_i + b*p_i^2
   dq_i = -(nu2 + b) * q_i + la2 * q_g + mu_i + b*q_i^2
   dq_g = -(mu_g + b) * q_g + mu_g + b*q_i*q_g
   dt = (x[1] >= pars[length(pars)])
   dx = c(dH,dp_i,dq_i,dq_g,dt)
   return(list(dx))
}

#' Function to simulate the approximate protracted speciation process
#'
#' Simulating the protracted speciation process according to the approxiomate
#' model of Lambert et al. 2014. This function differs from pbd_sim that 1) it
#' requires that the speciation-initiation rate is the same for good and
#' incipient species, and 2) that it does not simulate the exact protracted
#' speciation process, but an approximation made by the coalescent point
#' process.
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b (=
#' la_1 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to mu_1 (= mu_g in Etienne & Rosindell 2012) =
#' extinction rate of good species \cr \code{pars[3]} corresponds to la_1 (=
#' la_2 in Etienne & Rosindell 2012) = speciation completion rate \cr
#' \code{pars[4]} corresponds to mu_2 (= mu_i in ER2012) = extinction rate of
#' incipient species \cr When rates depend on time this time dependence should
#' be specified in pars1f and pars1 then becomes the parameters used in pars1f
#' \cr \cr
#' @param parsf Vector of functions how the rates depend on time, default
#' functions are constant functions of the parameters in pars1: \cr \cr
#' \code{parsf[1]} corresponds to time-dependence of b (= la_1 in Etienne &
#' Rosindell R2012) = speciation initiation rate \cr \code{parsf[2]}
#' corresponds to time-dependence of mu_1 (= mu_g in Etienne & Rosindell 2012)
#' = extinction rate of good species \cr \code{parsf[3]} corresponds to
#' tiem-dependence of la_1 (= la_2 in Etienne & Rosindell 2012) = speciation
#' completion rate \cr \code{parsf[4]} corresponds to time-dependence of mu_2
#' (= mu_i in ER2012) = extinction rate of incipient species \cr \cr
#' @param age Sets the crown age for the simulation
#' @param ntips Set the number of tips. If NULL (default) the number of tips will
#' be sampled
#' @param soc Determines whether the simulation should start at stem (1) or
#' crown (2) age
#' @param plotltt Sets whether the lineage-through-time plot should be plotted
#' (1) or not (0)
#' @param methode Sets which method should be used in the ode-solver. Default
#' is 'lsoda'. See package deSolve for details.
#' @return A set of branching times
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_sim}}
#' @keywords models
#' @examples
#'  pbd_sim_cpp(pars = c(0.2,1,0.2,0.1),age = 15)
#' @export pbd_sim_cpp
pbd_sim_cpp = function(pars,
                       parsf = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}),
                       age,
                       ntips = NULL,
                       soc = 2,
                       plotltt = 1,
                       methode = "lsoda")
{
  pars1 = c(parsf,pars)
  abstol = 1e-16
  reltol = 1e-10
  probs = c(1,1,0,0,0)
  y = deSolve::ode(probs,c(0,age),pbd_loglik_rhs_cpp,c(pars1,0),rtol = reltol,atol = abstol,method = methode)
  pT = 1 - y[2,2]
  if(is.null(ntips)) {
    nd = sum(stats::rgeom(soc,1 - pT))
  } else {
    nd <- ntips - soc
  }
  brts = rep(0,nd + 1)
  brts[1] = age
  i = 1
  while(i <= nd)
  {
     probs = c(1,1,0,0,0)
     y = deSolve::ode(probs,c(0,age),pbd_loglik_rhs_cpp,c(pars1,1 - stats::runif(1)*pT),rtol = reltol,atol = abstol,method = methode)
     brts[i + 1] = y[2,6]
     i = i + 1
  }
  brts = sort(brts, decreasing = T)
  if(plotltt == 1)
  {
      graphics::plot(c(-brts,0),c(soc:(length(brts) + soc - 1),length(brts) + soc - 1),log = 'y',type = 's',xlab = 'Time',ylab = 'Number of lineages')
  }
  return(brts)
}
