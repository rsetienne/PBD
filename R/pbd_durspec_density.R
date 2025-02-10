#' Probability density for duration of speciation under protracted birth-death
#' model of diversification
#'
#' pbd_durspec_density computes the probability density of the duration of
#' speciation under the protracted speciation model for a given set of
#' parameters
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b (=
#' la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr \code{pars[3]} corresponds to mu_2 (= mu_i in
#' ER2012) = extinction rate of incipient species \cr
#' @param tau The duration of speciation for which the density must be computed
#' @return The probability density
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_durspec_cumdensity}}\cr
#' \code{\link{pbd_durspec_mean}}\cr \code{\link{pbd_durspec_mode}}\cr
#' \code{\link{pbd_durspec_quantile}}\cr \code{\link{pbd_durspec_moment}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#'  pbd_durspec_density(pars = c(0.5,0.3,0.1), tau = 1)
#' @export pbd_durspec_density
pbd_durspec_density = function(pars,tau)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
if(la2 < Inf)
{
   phi = la2 - la3 + mu2
   D = sqrt((la2 + la3)^2 + 2*(la2 - la3) * mu2 + mu2^2)
   rho = 2*D^2 * exp(-D*tau) * (D + phi) / (D + phi + exp(-D*tau) * (D - phi))^2
} else {
   rho = Inf * (tau == 0)
   rho[is.nan(rho)] = 0
}
return(rho)
}
