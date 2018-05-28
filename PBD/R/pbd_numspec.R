pbd_pgeom <- function(pars,parsf = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), age, soc = 2, methode = "lsoda")
{
  checkpars(pars,age,soc)
  pars1 <- c(parsf,pars)
  abstol <- 1e-16
  reltol <- 1e-10
  probs <- c(1,1,0,0,0)
  y <- deSolve::ode(probs,c(0,age),pbd_loglik_rhs_cpp,c(pars1,0),rtol = reltol,atol = abstol,method = methode)
  pT <- y[2,2]
  return(pT)
}

checkpars <- function(pars,age,soc)
{
  if(min(c(pars,age)) < 0)
  {
    stop('Parameters and clade age must be non-negative real numbers.')
  }
  if(soc != 1 && soc != 2)
  {
    stop('soc must be 1 or 2.')
  }

  return(NULL)
}

checkquantile <- function(quantile)
{
  if(quantile < 0 || quantile > 1)
  {
    stop('Quantile must be a real number between 0 and 1.')
  }
  return(NULL)
}

#' @title Mean number of species under protracted birth-death model of
#' diversification
#' @description pbd_numspec_mean computes the mean number of species under the
#' protracted speciation model for a given set of parameters
#' @param pars Vector of parameters: \cr \cr
#' \code{pars[1]} corresponds to b (=
#' la_1 = la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to mu_1 (= mu_g in
#' ER2012) = extinction rate of good species \cr
#' \code{pars[3]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr
#' \code{pars[4]} corresponds
#' to mu_2 (= mu_i in ER2012) = extinction rate of incipient species \cr
#' @param age the stem or crown age (see soc)
#' @param soc specify whether it is the stem or the crown age
#' @return The expected duration of speciation
#' @author Rampal S. Etienne
#' @seealso
#' \code{\link{pbd_numspec_density}}\cr
#' \code{\link{pbd_numspec_cumdensity}}\cr
#' \code{\link{pbd_numspec_quantile}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#' pbd_numspec_mean(pars = c(0.3,0.1,0.5,0.1), age = 10, soc = 2)
#' @export pbd_numspec_mean
pbd_numspec_mean <- function(pars,parsf = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), age, soc = 2, methode = "lsoda")
{
  pT <- pbd_pgeom(pars = pars,parsf = parsf,age = age,soc = soc,methode = methode)
  numspec_mean <- (soc == 2) + as.numeric(soc * (1 - pT)/pT)

  #to be put in a test
  #endmc = 1000000
  #nd = rep(0,endmc)
  #for(mc in 1:endmc)
  #{
  #   nd[mc] = (soc == 2) + sum(stats::rgeom(soc,pT))
  #}
  #print(mean(nd))
  return(numspec_mean)
}



#' @title Number of species for a given quantile under protracted birth-death model of
#' diversification
#' @description pbd_numspec_quantile computes the number of species for a given quantile under the
#' protracted speciation model for a given set of parameters
#' @param pars Vector of parameters: \cr \cr
#' \code{pars[1]} corresponds to b (=
#' la_1 = la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to mu_1 (= mu_g in
#' ER2012) = extinction rate of good species \cr
#' \code{pars[3]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr
#' \code{pars[4]} corresponds
#' to mu_2 (= mu_i in ER2012) = extinction rate of incipient species \cr
#' @param age the stem or crown age (see soc)
#' @param soc specify whether it is the stem or the crown age
#' @return The expected duration of speciation
#' @author Rampal S. Etienne
#' @seealso
#' \code{\link{pbd_numspec_density}}\cr
#' \code{\link{pbd_numspec_cumdensity}}\cr
#' \code{\link{pbd_numspec_quantile}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#' pbd_numspec_quantile(pars = c(0.3,0.1,0.5,0.1), age = 10, soc = 2, quantile = 0.95)
#' @export pbd_numspec_quantile
pbd_numspec_quantile <- function(pars,parsf = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), age, soc = 2, methode = "lsoda", quantile)
{
  pT <- pbd_pgeom(pars = pars,parsf = parsf,age = age,soc = soc,methode = methode)
  numspec_quantile <- (soc == 2) + qnbinom(size = soc,prob = pT, p = quantile)
  return(numspec_quantile)
}



#' @title Median number of species under protracted birth-death model of
#' diversification
#' @description pbd_numspec_median computes the median number of species under the
#' protracted speciation model for a given set of parameters
#' @param pars Vector of parameters: \cr \cr
#' \code{pars[1]} corresponds to b (=
#' la_1 = la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to mu_1 (= mu_g in
#' ER2012) = extinction rate of good species \cr
#' \code{pars[3]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr
#' \code{pars[4]} corresponds
#' to mu_2 (= mu_i in ER2012) = extinction rate of incipient species \cr
#' @param age the stem or crown age (see soc)
#' @param soc specify whether it is the stem or the crown age
#' @return The expected duration of speciation
#' @author Rampal S. Etienne
#' @seealso
#' \code{\link{pbd_numspec_density}}\cr
#' \code{\link{pbd_numspec_cumdensity}}\cr
#' \code{\link{pbd_numspec_quantile}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#' pbd_numspec_median(pars = c(0.3,0.1,0.5,0.1), age = 10, soc = 2)
#' @export pbd_numspec_median
pbd_numspec_median <- function(pars,parsf = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), age, soc = 2, methode = "lsoda")
{
   return(pbd_numspec_quantile(pars = pars,parsf = parsf, age = age, soc = soc, methode = methode, quantile = 0.5))
}
