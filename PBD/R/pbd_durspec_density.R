<<<<<<< HEAD
#' @export

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
=======
#' @export

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
>>>>>>> 02bbac1e016a4e33fb6caad62db78bceaf1d65bc
