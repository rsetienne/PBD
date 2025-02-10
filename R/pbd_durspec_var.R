#' Variance in duration of speciation under protracted birth-death model of
#' diversification
#'
#' pbd_durspec_var computes the variance in the duration of speciation under
#' the protracted speciation model for a given set of parameters
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b (=
#' la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr \code{pars[3]} corresponds to mu_2 (= mu_i in
#' ER2012) = extinction rate of incipient species \cr
#' @return The variance in the duration of speciation
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_durspec_density}}\cr
#' \code{\link{pbd_durspec_cumdensity}}\cr \code{\link{pbd_durspec_mean}}\cr
#' \code{\link{pbd_durspec_mode}}\cr \code{\link{pbd_durspec_quantile}}\cr
#' \code{\link{pbd_durspec_moment}}
#' @keywords models
#' @examples
#'  pbd_durspec_var(pars = c(0.5,0.3,0.1))
#' @export pbd_durspec_var
pbd_durspec_var = function(pars)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
rho_var = pbd_durspec_moment(pars,2) - (pbd_durspec_mean(pars))^2
return(rho_var)
}
