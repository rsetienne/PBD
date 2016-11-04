<<<<<<< HEAD
#' @export

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
=======
#' @export

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
>>>>>>> 02bbac1e016a4e33fb6caad62db78bceaf1d65bc
