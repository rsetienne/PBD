#' mode of the duration of speciation under protracted birth-death model of
#' diversification
#'
#' pbd_durspec_mode computes the mode of the duration of speciation under the
#' protracted speciation model for a given set of parameters
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b (=
#' la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr \code{pars[3]} corresponds to mu_2 (= mu_i in
#' ER2012) = extinction rate of incipient species \cr
#' @return The expected duration of speciation
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_durspec_density}}\cr
#' \code{\link{pbd_durspec_cumdensity}}\cr \code{\link{pbd_durspec_mean}}\cr
#' \code{\link{pbd_durspec_quantile}}\cr \code{\link{pbd_durspec_moment}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#'  pbd_durspec_mode(pars = c(0.5,0.3,0.1))
#' @export pbd_durspec_mode
pbd_durspec_mode = function(pars)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
phi = la2 - la3 + mu2
D = sqrt((la2 + la3)^2 + 2*(la2 - la3) * mu2 + mu2^2)
tau_mode = max(0,1/D * log((D - phi)/(D + phi)))
return(tau_mode)
}
