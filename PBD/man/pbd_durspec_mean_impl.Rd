% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbd_durspec_mean.R
\name{pbd_durspec_mean_impl}
\alias{pbd_durspec_mean_impl}
\title{Actual calculation of the mean duration of speciation (equations 19 and 20
of reference article) assuming all inputs are correct}
\usage{
pbd_durspec_mean_impl(la2, la3, mu2)
}
\arguments{
\item{la2}{lambda_2, the speciation completion rate, in probability per time
unit}

\item{la3}{lambda_3, speciation initiation rate of incipient species, in
probability per time unit}

\item{mu2}{mu_2 extinction rate of the incipient species, in probability per
time unit}
}
\description{
Actual calculation of the mean duration of speciation (equations 19 and 20
of reference article) assuming all inputs are correct
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
Rampal S. Etienne
}
