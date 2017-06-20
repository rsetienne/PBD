#' Calculate the mean duration of speciation
#' @param pars
#' @export
pbd_durspec_mean = function(pars)
{
  la3 = pars[1]
  la2 = pars[2]
  mu2 = pars[3]
  if(la2 == Inf)
  {
    rho_mean = 0
  } else if(la2 == 0)
  {
    rho_mean = Inf
  } else if(la3 == 0)
  {
    rho_mean = 1/(la2 + mu2)
  } else if(mu2 == 0)
  {
    rho_mean = 1/la3 * log(1 + la3/la2)
  } else
  {
    D = sqrt((la2 + la3)^2 + 2*(la2 - la3) * mu2 + mu2^2)
    rho_mean = 2/(D - la2 + la3 - mu2) * log(2 / (1 + (la2 - la3 + mu2)/D))
  }
  return(rho_mean)
}

#' Calculate the mean duration of speciation
#' @export
pbd_mean_durspec = function(pars)
