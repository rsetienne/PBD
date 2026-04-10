#' @title Mean number of good species with n incipient species and of n orphaned
#' incipient species under protracted birth-death model of diversification
#' @description pbd_numspec_vec2 computes the mean number of good species with n
#' incipient daughter species and the number of cases with n orphaned incipient
#' siblings under the protracted speciation model for a given set of parameters
#' @param pars Vector of parameters: \cr \cr
#' \code{pars[1]} corresponds to b_G (=
#' la_1 in Etienne & Rosindell R2012) = speciation initiation rate of good
#' species\cr
#' \code{pars[2]} corresponds to b_I (=
#' la_3 in Etienne & Rosindell R2012) = speciation initiation rate of incipient
#' species\cr
#' \code{pars[3]} corresponds to mu_G = extinction rate of good species \cr
#' \code{pars[4]} corresponds to mu_I = extinction rate of incipient species \cr
#' \code{pars[5]} corresponds to la (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr
#' @param age the stem age
#' @param initvec initial vector. Default is 1 good species with 0 incipient
#' species.
#' @return The expected number of good species with n incipient specis and the number of
#' caes with n orphaned incipient siblings
#' @author Rampal S. Etienne
#' @seealso
#' \code{\link{pbd_numspec_mean}}\cr
#' \code{\link{pbd_numspec_quantile}}\cr
#' \code{\link{pbd_numspec_median}}\cr
#' \code{\link{pbd_numspec_vec}}
#' @keywords models
#' @examples
#' pbd_numspec_vec2(pars = c(0.3,0.3,0.1,0.1,0.2), age = 10)
#' @export pbd_numspec_vec2
pbd_numspec_vec2 <- function(pars = c(0.3,0.3,0.1,0.1,100),
                             age = 10,
                             initvec = c(1,rep(0,20)),
                             abstol = 1E-10,
                             reltol = 1E-10,
                             methode = 'odeint::runge_kutta_cash_karp54') {
  if(length(initvec) %% 2 == 0) stop('The length of the initial vector should be odd')
  distr <- pbd_integrate_odeint(initvec,
                                c(-age,0),
                                pars,
                                abstol,
                                reltol,
                                methode)
  lx <- length(distr)
  good <- distr[1:((lx + 1)/2)]
  orphaned <- distr[((lx + 1)/2 + 1):lx]
  numspec_tot <- sum(good * c(1:((lx + 1)/2))) +
                 sum(orphaned * c(1:((lx + 1)/2 - 1)))
  numspec_rep <- sum(distr)
  return(list(good = good,
              orphaned = orphaned,
              numspec_tot = numspec_tot,
              numspec_rep = numspec_rep))
}
