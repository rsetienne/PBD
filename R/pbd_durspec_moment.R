#' Moments of duration of speciation under protracted birth-death model of
#' diversification
#'
#' pbd_durspec_moment computes the moments of the duration of speciation under
#' the protracted speciation model for a given set of parameters
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b (=
#' la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr \code{pars[3]} corresponds to mu_2 (= mu_i in
#' ER2012) = extinction rate of incipient species \cr
#' @param order order of the moment to compute (1 is first moment, giving the
#' mean)
#' @return The moment of the duration of speciation
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_durspec_density}}\cr
#' \code{\link{pbd_durspec_cumdensity}}\cr \code{\link{pbd_durspec_mean}}\cr
#' \code{\link{pbd_durspec_mode}}\cr \code{\link{pbd_durspec_quantile}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#'  pbd_durspec_moment(pars = c(0.5,0.3,0.1),2)
#' @export pbd_durspec_moment
pbd_durspec_moment = function(pars,order)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
if(la2 < Inf)
{
   rho_moment = stats::integrate(function(x) {x^order * pbd_durspec_density(pars,x) },lower = 0, upper = Inf, abs.tol = 1e-10)
} else {
   if(order == 0)
   {
       rho_moment = 1
   } else {
       rho_moment = 0
   }
}
return(rho_moment$value)
}
