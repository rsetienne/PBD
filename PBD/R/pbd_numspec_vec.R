#' @title Mean number of good and incipient species under protracted birth-death
#' model of diversification
#' @description pbd_numspec_vec computes the mean number of good and incipient
#' species under the protracted speciation model for a given set of parameters
#' @param pars Vector of parameters: \cr \cr
#' \code{pars[1]} corresponds to b_G (=
#' la_1 in Etienne & Rosindell R2012) = speciation initiation rate of good
#' species\cr
#' \code{pars[2]} corresponds to b_G (=
#' la_1 in Etienne & Rosindell R2012) = speciation initiation rate of incipient
#' pecies\cr
#' \code{pars[3]} corresponds to mu_G = extinction rate of good species \cr
#' \code{pars[4]} corresponds to mu_I = extinction rate of incipient species \cr
#' \code{pars[5]} corresponds to la (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr
#' @param age the stem age
#' @param initvec initial vector of good and incipient species. Default is 1 good
#' species.
#' @details The result follows from a simple 2D system of differential equations
#' for G and I.
#' @return The expected number of representative species
#' @author Rampal S. Etienne
#' @seealso
#' \code{\link{pbd_numspec_mean}}\cr
#' \code{\link{pbd_numspec_quantile}}\cr
#' \code{\link{pbd_numspec_median}}
#' @keywords models
#' @examples
#' pbd_numspec_vec(pars = c(0.3,0.3,0.1,0.1,0.2), age = 10)
#' @export pbd_numspec_vec
pbd_numspec_vec <- function(pars, age, initvec = c(1,0)) {
  b_G <- pars[1]
  b_I <- pars[2]
  mu_G <- pars[3]
  mu_I <- pars[4]
  lambda <- pars[5]
  M <- matrix(c(-mu_G, lambda, b_G, b_I - lambda - mu_I),2,2, byrow = TRUE)
  return(Matrix::expm(M * age) %*% initvec)
}
