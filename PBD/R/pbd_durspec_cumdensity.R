<<<<<<< HEAD
#' @export

pbd_durspec_cumdensity = function(pars,tau)
{
   stats::integrate(function(x) pbd_durspec_density(pars,x),lower = 0,upper = tau, abs.tol = 1e-10)$value
}
=======
#' @export

pbd_durspec_cumdensity = function(pars,tau)
{
   stats::integrate(function(x) pbd_durspec_density(pars,x),lower = 0,upper = tau, abs.tol = 1e-10)$value
}
>>>>>>> 02bbac1e016a4e33fb6caad62db78bceaf1d65bc
