% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbd_durspec_mean.R
\name{pbd_mean_durspecs}
\alias{pbd_mean_durspecs}
\title{Calculate the mean durations of speciation (equations 19 and 20 of reference
article)}
\usage{
pbd_mean_durspecs(eris, scrs, siris)
}
\arguments{
\item{eris}{one or more extinction rates of the incipient species, or mu_2
in article, in probability per time unit. These values will be recycled if
needed}

\item{scrs}{one or more speciation completion rates, or lambda_2 in article,
in probability per time unit. These values will be recycled if needed}

\item{siris}{one or more speciation initiation rates of incipient species,
or lambda_3 in article, in probability per time unit. These values will be
recycled if needed}
}
\value{
the means durations of speciation, in time units. Puts an NA at each
invalid combination of inputs
}
\description{
Calculate the mean durations of speciation (equations 19 and 20 of reference
article)
}
\examples{

  eris <- c(0.1, 0.2) # extinction rates of incipient species
  scrs <- c(0.2, 0.3)  # speciation completion rates
  siris <- c(0.3, 0.4) # speciation initiation rates of incipient species
  mean_durspecs <- pbd_mean_durspecs(eris, scrs, siris)
  expected_mean_durspecs <- c(2.829762, 1.865386)
  testthat::expect_equal(mean_durspecs, expected_mean_durspecs,
    tolerance = 0.000001)

}
\references{
Etienne, Rampal S., and James Rosindell. "Prolonging the past
counteracts the pull of the present: protracted speciation can explain
observed slowdowns in diversification." Systematic Biology 61.2 (2012):
204-213.
}
\seealso{
pbd_mean_durspec
}
\author{
Richel J.C. Bilderbeek
}
