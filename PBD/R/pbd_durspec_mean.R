#' Calculate the mean duration of speciation
#' (equations 19 and 20 of reference article),
#' unwraps its arguments without checking for input
#' @param pars numeric vector with three elements:
#'   [1] speciation initiation rate of incipient species (lambda_3),
#'   [2] speciation completion rate (lamba_2) and
#'   [3] extinction rate of the incipient species (mu_2)
#' @return the means duration of speciation, in time units
#' @examples
#'   eri <- 0.1 # extinction rate of incipient species
#'   scr <- 0.2 # speciation completion rate
#'   siri <- 0.3 # speciation initiation rate of incipient species
#'   mean_durspec <- pbd_durspec_mean(c(siri, scr, eri))
#'   expected_mean_durspec <- 2.829762
#'   testthat::expect_equal(mean_durspec, expected_mean_durspec,
#'     tolerance = 0.000001)
#' @author Rampal S. Etienne, Richel J.C. Bilderbeek
#' @references Etienne, Rampal S., and James Rosindell. "Prolonging the past
#'   counteracts the pull of the present: protracted speciation can explain
#'   observed slowdowns in diversification." Systematic
#'   Biology 61.2 (2012): 204-213.
#' @seealso pbd_mean_durspec
#' @note future releases may check if user input is valid
#' @export
pbd_durspec_mean = function(pars)
{
  # Do not check 'pars' being valid, as this is the classic behavior
  pbd_durspec_mean_impl(la2 = pars[2], la3 = pars[1], mu2 = pars[3])
}

#' Calculate the mean durations of speciation
#' (equations 19 and 20 of reference article)
#' @param scrs one or more speciation completion rates,
#'   or lambda_2 in article,
#'   in probability per time unit.
#'   These values will be recycled if needed
#' @param siris one or more speciation initiation rates of incipient species, or
#'   lambda_3 in article,
#'   in probability per time unit.
#'   These values will be recycled if needed
#' @param eris one or more extinction rates of the incipient species,
#'   or mu_2 in article,
#'   in probability per time unit.
#'   These values will be recycled if needed
#' @return the means durations of speciation, in time units. Puts an NA
#'   at each invalid combination of inputs
#' @examples
#'   eris <- c(0.1, 0.2) # extinction rates of incipient species
#'   scrs <- c(0.2, 0.3)  # speciation completion rates
#'   siris <- c(0.3, 0.4) # speciation initiation rates of incipient species
#'   mean_durspecs <- pbd_mean_durspecs(eris, scrs, siris)
#'   expected_mean_durspecs <- c(2.829762, 1.865386)
#'   testthat::expect_equal(mean_durspecs, expected_mean_durspecs,
#'     tolerance = 0.000001)
#' @author Richel J.C. Bilderbeek
#' @references Etienne, Rampal S., and James Rosindell. "Prolonging the past
#'   counteracts the pull of the present: protracted speciation can explain
#'   observed slowdowns in diversification." Systematic
#'   Biology 61.2 (2012): 204-213.
#' @seealso pbd_mean_durspec
#' @export
pbd_mean_durspecs = function(eris, scrs, siris) {

  # Find invalid indices
  invalid <- eris < 0.0 | is.na(eris) |
    scrs < 0.0 | is.na(scrs) |
    siris < 0.0 | is.na(siris)

  # Correct invalid rates to valid ones
  correct <- function(x) {
    x[ is.na(x) | x < 0.0] <- 0.0
    x
  }

  # Get durations for corrected rates
  v <- mapply(pbd_mean_durspec, correct(eris), correct(scrs), correct(siris))

  # Let invalid rates become NAs
  v[ invalid ] <- NA
  v
}

#' Calculate the mean duration of speciation
#' (equations 19 and 20 of reference article), non-vectorized
#' @param scr one single speciation completion rate,
#'   or lambda_2 in article,
#'   in probability per time unit
#' @param siri one single speciation initiation rate of incipient species, or
#'   lambda_3 in article,
#'   in probability per time unit
#' @param eri one single extinction rate of the incipient species,
#'   or mu_2 in article,
#'   in probability per time unit
#' @return the means duration of speciation, in time units
#' @examples
#'   eri <- 0.1 # extinction rate of incipient species
#'   scr <- 0.2 # speciation completion rate
#'   siri <- 0.3 # speciation initiation rate of incipient species
#'   mean_durspec <- pbd_mean_durspec(eri, scr, siri)
#'   expected_mean_durspec <- 2.829762
#'   testthat::expect_equal(mean_durspec, expected_mean_durspec,
#'     tolerance = 0.000001)
#' @author Richel J.C. Bilderbeek
#' @references Etienne, Rampal S., and James Rosindell. "Prolonging the past
#'   counteracts the pull of the present: protracted speciation can explain
#'   observed slowdowns in diversification." Systematic
#'   Biology 61.2 (2012): 204-213.
#' @export
pbd_mean_durspec = function(eri, scr, siri) {
  if (is.na(eri) || eri < 0.0) {
    stop("extinction rate of incipient species must be zero or positive")
  }
  if (is.na(scr) || scr < 0.0) {
    stop("speciation completion rate must be zero or positive")
  }
  if (is.na(siri) || siri < 0.0) {
    stop("speciation initiation rate of incipient species must ",
      "be zero or positive")
  }
  pbd_durspec_mean_impl(la2 = scr, la3 = siri, mu2 = eri)
}


#' Actual calculation of the mean duration of speciation
#' (equations 19 and 20 of reference article)
#' assuming all inputs are correct
#' @param la2 lambda_2, the speciation completion rate,
#'   in probability per time unit
#' @param la3 lambda_3, speciation initiation rate of incipient species,
#'   in probability per time unit
#' @param mu2 mu_2 extinction rate of the incipient species,
#'   in probability per time unit
#' @author Rampal S. Etienne
#' @references Etienne, Rampal S., and James Rosindell. "Prolonging the past
#'   counteracts the pull of the present: protracted speciation can explain
#'   observed slowdowns in diversification." Systematic
#'   Biology 61.2 (2012): 204-213.
#' @seealso pbd_mean_durspec
pbd_durspec_mean_impl = function(la2, la3, mu2)
{
  if(la2 == Inf) {
    rho_mean <- 0.0
  } else if(la2 == 0) {
    rho_mean <- Inf
  } else if(la3 == 0) {
    rho_mean <- 1.0 / (la2 + mu2)
  } else if(mu2 == 0) {
    rho_mean <- 1/la3 * log(1 + la3/la2)
  } else {
    D <- sqrt((la2 + la3)^2 + 2*(la2 - la3) * mu2 + mu2^2)
    rho_mean <- 2/(D - la2 + la3 - mu2) * log(2 / (1 + (la2 - la3 + mu2)/D))
  }
  rho_mean
}
