

#' Protracted Birth-Death Model of Diversification
#' 
#' This package computes the (maximum) likelihood of the protracted speciation
#' model for a given set of branching times This package is a likelihood-based
#' statistical package to estimate parameters under the protracted speciation
#' model.\cr \cr First version: 0.8\cr New in version 0.9\cr - Bug fix for stem
#' age\cr New in version 0.91\cr - Reports loglik = -Inf on an error in the
#' deSolve package (function ode)\cr New in version 0.92\cr - Correcting order
#' of parameters of pbd_sim\cr New in version 0.93\cr - pbd_sim produces a
#' tree, a matrix containing all events in the simulation, and a tree with one
#' sample per species.\cr New in version 1.0\cr - Conditioning is also possible
#' on a range of values of the number of species.\cr New in version 1.1\cr -
#' Simulation of the protracted speciation tree has more features. \cr New in
#' version 1.2\cr - Optimization can make use of subplex (default) and simplex
#' (older versions).\cr New in version 1.3\cr - Contains a function to carry
#' out a bootstrap likelihood ratio test.\cr - Vignette and test added.\cr -
#' Reports an error if exteq = TRUE and initparsopt contains 4 parameters.\cr -
#' Option to limit a simulation to a certain maximum number of species; if
#' exceeded, the simulation is ignored.\cr New in version 1.4:\cr - Includes
#' all special cases in pbd_durspec_mean\cr - Fixes a bug in conditioning on a
#' range of values of the number of species\cr
#' 
#' pbd_loglik computes the likelihood of the protracted birth-death model of
#' diversification, given a set of parameters and a data set of phylogenetic
#' branching times.
#' 
#' pbd_ML finds the parameters that maximizes the likelihood computed by
#' pbd_loglik.
#' 
#' pbd_bootstrap performs a maximum likelihood analysis and simulates with the
#' maximum likelihood parameters. The ML parameters of the simulated data sets
#' are then estimated, providing an uncertainty distribution for the original
#' ML estimate on the original data.
#' 
#' @name PBD-package
#' @aliases PBD-package PBD
#' @docType package
#' @author Rampal S. Etienne Maintainer: Rampal S. Etienne
#' <r.s.etienne@@rug.nl>
#' @seealso \code{DDD}
#' @references - Etienne, R.S. & J. Rosindell 2012. Systematic Biology 61:
#' 204-213.\cr - Lambert, A., H. Morlon & R.S. Etienne 2014. Journal of
#' Mathmematical Biology 70: 367-397. doi:10.1007/s00285-014-0767-x\cr -
#' Etienne, R.S., H. Morlon & A. Lambert 2014. Evolution 68: 2430-2440,
#' doi:10.1111/evo.12433\cr. - Etienne, R.S., A.L. Pigot & A.B. Phillimore
#' 2016. Methods in Ecology & Evolution 7: 1092-1099, doi:
#' 10.1111/2041-210X.12565 \cr
#' @keywords package
NULL



