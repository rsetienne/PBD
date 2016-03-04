\name{pbd_ML}
\alias{pbd_ML}
\title{Maximization of loglikelihood under protracted birth-death model of diversification}
\description{Likelihood maximization for protracted birth-death model of diversification}
\usage{
pbd_ML(
   brts,
   initparsopt = c(0.2,0.1,1),
   idparsopt = 1:length(initparsopt),
   idparsfix = NULL,
   parsfix = NULL,
   exteq = 1,
   parsfunc = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},
function(t,pars) {pars[3]}, function(t,pars) {pars[4]}),
   missnumspec = 0,
   cond = 1,
   btorph = 1,
   soc = 2,
   methode = "lsoda",
   n_low = 0,
   n_up = 0,
   tol = c(1E-6, 1E-6, 1E-6),
   maxiter = 1000 * round((1.25)^length(idparsopt)),
   optimmethod = 'subplex'
)}
\arguments{
  \item{brts}{A set of branching times of a phylogeny, all positive }
  \item{initparsopt}{The initial values of the parameters that must be optimized}
  \item{idparsopt}{
   The ids of the parameters that must be optimized, e.g. 1:4 for all parameters.
   The ids are defined as follows: \cr
   id == 1 corresponds to b (speciation-initiation rate) \cr
   id == 2 corresponds to mu_1 (extinction rate of good species) \cr
   id == 3 corresponds to la_1 (speciation-completion rate) \cr
   id == 4 corresponds to mu_2 (extinction rate of incipient species) \cr
   }
  \item{idparsfix}{
The ids of the parameters that should not be optimized, e.g. c(2,4) if mu_1 and mu_2
should not be optimized, but only b and la_1. In that case idparsopt must be c(1,3).
  }
  \item{parsfix}{ The values of the parameters that should not be optimized }
  \item{exteq}{Sets whether incipient species have the same (1) or different (0) extinction rate as good species.
If exteq = 0, then idparsfix and idparsopt should together have all parameters 1:4
}
  \item{parsfunc}{ Specifies functions how the rates depend on time, default functions are constant functions }
  \item{missnumspec}{
The number of species that are in the clade but missing in the phylogeny }
  \item{cond}{Conditioning: \cr
    cond == 0 : conditioning on stem or crown age \cr
    cond == 1 : conditioning on stem or crown age and non-extinction of the phylogeny \cr
    cond == 2 : conditioning on stem or crown age and number of extant taxa \cr
}
  \item{btorph}{Sets whether the likelihood is for the branching times (0) or the phylogeny (1)}
  \item{soc}{Sets whether the first element of the branching times is the stem (1) or the crown (2) age}
  \item{methode}{Sets which method should be used in the ode-solver. Default is 'lsoda'.
See package deSolve for details.}
  \item{n_low}{Sets the lower bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)}
  \item{n_up}{Sets the upper bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)}  
  \item{tol}{Sets the tolerances in the optimization. Consists of: \cr
   reltolx = relative tolerance of parameter values in optimization \cr
   reltolf = relative tolerance of function value in optimization \cr
   abstolx = absolute tolerance of parameter values in optimization
}
  \item{maxiter}{Sets the maximum number of iterations in the optimization}
  \item{optimmethod}{ Method used in optimization of the likelihood. Current default is 'subplex'. Alternative is 'simplex' (default of previous versions)}  
}
\value{
\item{b}{ gives the maximum likelihood estimate of b}
 \item{mu_1}{ gives the maximum likelihood estimate of mu_1}
 \item{la_1}{ gives the maximum likelihood estimate of la_1}
 \item{mu_2}{ gives the maximum likelihood estimate of mu_2}
 \item{loglik}{ gives the maximum loglikelihood}
 \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom}
 \item{conv}{ gives a message on convergence of optimization; conv = 0 means convergence}
}
\author{
Rampal S. Etienne
}
\seealso{
 \code{\link{pbd_loglik}}
}
\examples{
pbd_ML(1:10,initparsopt = c(0.2,0.01,0.3), exteq = 1)
}
\keyword{ models }
