<<<<<<< HEAD
#' @export

pbd_durspec_var = function(pars)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
rho_var = pbd_durspec_moment(pars,2) - (pbd_durspec_mean(pars))^2
return(rho_var)
}
=======
#' @export

pbd_durspec_var = function(pars)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
rho_var = pbd_durspec_moment(pars,2) - (pbd_durspec_mean(pars))^2
return(rho_var)
}
>>>>>>> 02bbac1e016a4e33fb6caad62db78bceaf1d65bc
