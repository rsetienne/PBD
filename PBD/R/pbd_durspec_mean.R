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
#' (equations 19 and 20 of reference article)
#' @param scr speciation completion rate, in probability per time unit
#' @param siri speciation initiation rate of incipient species, in
#'   probability per time unit
#' @param eri extinction rate of the incipient species,
#'   in probability per time unit
#' @return the means duration of speciation, in time units
#' @examples
#'   eri <- 0.1 # incipient species extinction rate
#'   scr <- 0.2 # speciation completion rate
#'   siri <- 0.3 # speciation initiation rate of incipient species
#'   mean_durspec <- pbd_mean_durspec(eri, scr, siri)
#'   expected_mean_durspec <- 2.829762
#'   testthat::expect_equal(mean_durspec, expected_mean_durspec,
#'     tolerance = 0.000001)
#' @references Etienne, Rampal S., and James Rosindell. "Prolonging the past
#'   counteracts the pull of the present: protracted speciation can explain
#'   observed slowdowns in diversification." Systematic
#'   Biology 61.2 (2012): 204-213.
#' @seealso pbd_durspec_mean
#' @export
pbd_mean_durspec = function(eri, scr, siri) {
  if (is.na(scr) || scr < 0.0) {
    stop("speciation completion rate must be zero or positive")
  }
  if (is.na(siri) || siri < 0.0) {
    stop("speciation initiation rate of incipient species must",
      "be zero or positive")
  }
  if (is.na(eri) || eri < 0.0) {
    stop("extinction rate of incipient species must be zero or positive")
  }
  pbd_durspec_mean(c(siri, scr, eri))
}
