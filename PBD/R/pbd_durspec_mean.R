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
#' @param lambda_2 speciation completion rate, in probability per time unit
#' @param lambda_3 speciation initiation rate of incipient species, in
#'   probability per time unit
#' @param mu_2 incipient species extinction rate, in probability per time unit
#' @return the means duration of speciation, in time units
#' @examples
#'   lambda_2 <- 0.1 # speciation completion rate
#'   lambda_3 <- 0.2 # speciation initiation rate of incipient species
#'   mu_2 <- 0.3 # incipient species extinction rate
#'   mean_durspec <- pbd_mean_durspec(lambda_2, lambda_3, mu_2)
#'   expected_mean_durspec <- 3.242955
#'   tolerance <- 0.0000001
#'   if (abs(mean_durspec - expected_mean_durspec) > tolerance) {
#'     stop("Duration of speciation must match expectation")
#'   }
#' @references Etienne, Rampal S., and James Rosindell. "Prolonging the past
#'   counteracts the pull of the present: protracted speciation can explain
#'   observed slowdowns in diversification." Systematic
#'   Biology 61.2 (2012): 204-213.
#' @seealso pbd_durspec_mean
#' @export
pbd_mean_durspec = function(lambda_2, lambda_3, mu_2) {
  if (lambda_2 < 0.0) {
    stop("lambda_2 (speciation completion rate) cannot be negative")
  }
  if (lambda_3 < 0.0) {
    stop("lambda_3 (speciation initiation rate of incipient species) cannot",
      "be negative")
  }
  if (mu_2 < 0.0) {
    stop("mu_2 (incipient species extinction rate) cannot be negative")
  }
  pbd_durspec_mean(c(lambda_3, lambda_2, mu_2))
}
