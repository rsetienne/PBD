#' Quantiles of duration of speciation under protracted birth-death model of
#' diversification
#'
#' pbd_durspec_quantile computes a quantile of the duration of speciation under
#' the protracted speciation model for a given set of parameters
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b (=
#' la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars[2]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) =
#' speciation completion rate \cr \code{pars[3]} corresponds to mu_2 (= mu_i in
#' ER2012) = extinction rate of incipient species \cr
#' @param p Quantile (e.g. p = 0.5 gives the median)
#' @return The quantil of the duration of speciation
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_durspec_density}}\cr
#' \code{\link{pbd_durspec_cumdensity}}\cr \code{\link{pbd_durspec_mean}}\cr
#' \code{\link{pbd_durspec_mode}}\cr \code{\link{pbd_durspec_moment}}\cr
#' \code{\link{pbd_durspec_var}}
#' @keywords models
#' @examples
#'  pbd_durspec_quantile(pars = c(0.5,0.3,0.1),0.5)
#' @export pbd_durspec_quantile
pbd_durspec_quantile = function(pars,p)
{
   expdurspec = pbd_durspec_mean(pars)
   if(expdurspec < 1E-7)
   {
      q = 0
   } else {
      found = 0
      uptau = 100 * expdurspec
      while(found == 0)
      {
          if(pbd_durspec_cumdensity(pars,uptau) > p)
          {
              found = 1
          } else {
              uptau = 10*uptau
          }
      }
      q = stats::uniroot(function(x) pbd_durspec_cumdensity(pars,x) - p,c(0,uptau))$root
   }
   return(q)
}
