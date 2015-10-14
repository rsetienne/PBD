\name{PBD-package}
\alias{PBD-package}
\alias{PBD}
\docType{package}
\title{ Protracted birth-death model of diversification }
\description{ This package computes the (maximum) likelihood of the protracted speciation model for a given set of branching times
This package is a likelihood-based statistical package to estimate parameters
under the protracted speciation model.\cr \cr
   First version: 0.8\cr
   New in version 0.9\cr
   - Bug fix for stem age\cr
   New in version 0.91\cr
   - Reports loglik = -Inf on an error in the deSolve package (function ode)\cr
   New in version 0.92\cr
   - Correcting order of parameters of pbd_sim\cr
   New in version 0.93\cr
   - pbd_sim produces a tree, a matrix containing all events in the simulation, and a tree with one sample per species.\cr
   New in version 1.0\cr
   - Conditioning is also possible on a range of values of the number of species.\cr
   New in version 1.1\cr
   - Simulation of the protracted speciation tree has more features. \cr
}
\details{
\tabular{ll}{
Package: \tab PBD\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-9-16\cr
License: \tab GPL-2\cr
}
pbd_loglik computes the likelihood of the protracted birth-death model of
diversification, given a set of parameters and a data set of phylogenetic branching times.

pbd_ML finds the parameters that maximizes the likelihood computed by pbd_loglik.
}
\author{
Rampal S. Etienne
Maintainer: Rampal S. Etienne <r.s.etienne@rug.nl>
}
\references{
- Etienne, R.S. & J. Rosindell 2012. Systematic Biology 61: 204-213.\cr
- Lambert, A., H. Morlon & R.S. Etienne 2014. Journal of Mathmematical Biology 70: 367-397. doi:10.1007/s00285-014-0767-x\cr
- Etienne, R.S., H. Morlon & A. Lambert 2014. Evolution 68: 2430-2440 doi:10.1111/evo.12433\cr.
}
\keyword{ package }
\seealso{ \code{DDD} }
\examples{
pbd_ML(1:10)
}
