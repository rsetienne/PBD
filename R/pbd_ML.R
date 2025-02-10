#' Maximization of loglikelihood under protracted birth-death model of
#' diversification
#'
#' Likelihood maximization for protracted birth-death model of diversification
#'
#'
#' @param brts A set of branching times of a phylogeny, all positive
#' @param initparsopt The initial values of the parameters that must be
#' optimized
#' @param idparsopt The ids of the parameters that must be optimized, e.g. 1:4
#' for all parameters.  The ids are defined as follows: \cr id == 1 corresponds
#' to b (speciation-initiation rate) \cr id == 2 corresponds to mu_1
#' (extinction rate of good species) \cr id == 3 corresponds to la_1
#' (speciation-completion rate) \cr id == 4 corresponds to mu_2 (extinction
#' rate of incipient species) \cr
#' @param idparsfix The ids of the parameters that should not be optimized,
#' e.g. c(2,4) if mu_1 and mu_2 should not be optimized, but only b and la_1.
#' In that case idparsopt must be c(1,3).
#' @param parsfix The values of the parameters that should not be optimized
#' @param exteq Sets whether incipient species have the same (1) or different
#' (0) extinction rate as good species. If exteq = 0, then idparsfix and
#' idparsopt should together have all parameters 1:4
#' @param parsfunc Specifies functions how the rates depend on time, default
#' functions are constant functions
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param cond Conditioning: \cr cond == 0 : conditioning on stem or crown age
#' \cr cond == 1 : conditioning on stem or crown age and non-extinction of the
#' phylogeny \cr cond == 2 : conditioning on stem or crown age and number of
#' extant taxa \cr
#' @param btorph Sets whether the likelihood is for the branching times (0) or
#' the phylogeny (1)
#' @param soc Sets whether the first element of the branching times is the stem
#' (1) or the crown (2) age
#' @param methode Sets which method should be used in the ode-solver. Default
#' is 'lsoda'. See package deSolve for details.
#' @param n_low Sets the lower bound of the number of species on which
#' conditioning should be done when cond = 2. Set this to 0 when conditioning
#' should be done on precisely the number of species (default)
#' @param n_up Sets the upper bound of the number of species on which
#' conditioning should be done when cond = 2. Set this to 0 when conditioning
#' should be done on precisely the number of species (default)
#' @param tol Sets the tolerances in the optimization. Consists of: \cr reltolx
#' = relative tolerance of parameter values in optimization \cr reltolf =
#' relative tolerance of function value in optimization \cr abstolx = absolute
#' tolerance of parameter values in optimization
#' @param maxiter Sets the maximum number of iterations in the optimization
#' @param optimmethod Method used in optimization of the likelihood. Current
#' default is 'subplex'. Alternative is 'simplex' (default of previous
#' versions)
#' @param num_cycles Number of cycles of the optimization (default is 1).
#' @param verbose if TRUE, explanatory text will be shown
#' @return A data frame with the following components:\cr
#' \item{b}{ gives the maximum likelihood estimate of b}
#' \item{mu_1}{ gives the maximum likelihood estimate of mu_1}
#' \item{la_1}{ gives the maximum likelihood estimate of la_1}
#' \item{mu_2}{ gives the maximum likelihood estimate of mu_2}
#' \item{loglik}{ gives the maximum loglikelihood}
#' \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom}
#' \item{conv}{ gives a message on convergence of optimization;
#' conv = 0 means convergence}
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_loglik}}
#' @keywords models
#' @examples
#'
#' pbd_ML(1:10,initparsopt = c(4.640321,4.366528,0.030521), exteq = 1)
#'
#' @export pbd_ML
pbd_ML = function(
  brts,
  initparsopt = c(0.2,0.1,1),
  idparsopt = 1:length(initparsopt),
  idparsfix = NULL,
  parsfix = NULL,
  exteq = 1,
  parsfunc = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}),
  missnumspec = 0,
  cond = 1,
  btorph = 1,
  soc = 2,
  methode = "lsoda",
  n_low = 0,
  n_up = 0,
  tol = c(1E-6, 1E-6, 1E-6),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  optimmethod = 'subplex',
  num_cycles = 1,
  verbose = TRUE)
{
  #options(warn=-1)
  brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
  if(is.numeric(brts) == FALSE)
  {
    cat("The branching times should be numeric.\n")
    out2 = data.frame(b = -1, mu_1 = -1, lambda_1 = -1,mu_2 = -1, loglik = -1, df = -1, conv = -1)
  } else {
    if(exteq == 1){ idexteq = 4 } else { idexteq = NULL }
    idpars = sort(c(idparsopt,idparsfix,idexteq))
    if((prod(idpars == (1:4)) != 1) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
    {
      cat("The arguments should be coherent.\n")
      out2 = data.frame(b = -1, mu_1 = -1, lambda_1 = -1,mu_2 = -1, loglik = -1, df = -1, conv = -1)
    } else {
      namepars = c("b", "mu_1", "lambda_1", "mu_2")
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      if (verbose) { cat("You are optimizing",optstr,"\n") }
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      if (verbose) { cat("You are fixing",fixstr,"\n") }
      if(exteq == 1) { fixstr = "exactly" } else { fixstr = "not" }
      if (verbose) { cat("Extinction rate of incipient species is",fixstr,"the same as for good species.\n") }
      trparsopt = initparsopt/(1 + initparsopt)
      trparsfix = parsfix/(1 + parsfix)
      trparsfix[parsfix == Inf] = 1
      pars2 = c(cond,btorph,soc,0,methode,n_low,n_up)
      utils::flush.console()
      initloglik = pbd_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,exteq = exteq,parsfunc = parsfunc,pars2 = pars2,brts = brts,missnumspec = missnumspec)
      if (verbose) { cat("The likelihood for the initial parameter values is",initloglik,"\n") }
      utils::flush.console()
      if(initloglik == -Inf)
      {
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 = data.frame(b = -1, mu_1 = -1, lambda_1 = -1,mu_2 = -1, loglik = -1, df = -1, conv = -1)
      } else {
        if (verbose) { cat("Optimizing the likelihood - this may take a while.","\n") }
        utils::flush.console()
        optimpars = c(tol,maxiter)
        out = DDD::optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = pbd_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,exteq = exteq, parsfunc = parsfunc, pars2 = pars2,brts = brts, missnumspec = missnumspec, num_cycles = num_cycles)
        if(out$conv > 0)
        {
          cat("Optimization has not converged. Try again with different initial values.\n")
          out2 = data.frame(b = -1, mu_1 = -1, lambda_1 = -1,mu_2 = -1, loglik = -1, df = -1, conv = -1)
        } else {
          MLtrpars = as.numeric(unlist(out$par))
          MLpars = MLtrpars/(1-MLtrpars)
          MLpars1 = rep(0,4)
          MLpars1[idparsopt] = MLpars
          if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
          if(exteq == 1) { MLpars1[4] = MLpars[2] }
          if(MLpars1[3] > 10^7){MLpars1[3] = Inf}
          ML = as.numeric(unlist(out$fvalues))
          out2 = data.frame(b = MLpars1[1],mu_1 = MLpars1[2],lambda_1 = MLpars1[3], mu_2 = MLpars1[4], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
          if(verbose)
          {
             s1 = sprintf('Maximum likelihood parameter estimates: b: %f, mu_1: %f, lambda_1: %f, mu_2: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4])
             s2 = sprintf('Maximum loglikelihood: %f',ML)
             s3 = sprintf('The expected duration of speciation for these parameters is: %f',pbd_durspec_mean(c(MLpars1[1],MLpars1[3],MLpars1[4])))
             s4 = sprintf('The median duration of speciation for these parameters is: %f',pbd_durspec_quantile(c(MLpars1[1],MLpars1[3],MLpars1[4]),0.5))
             cat("\n",s1,"\n",s2,"\n",s3,"\n",s4,"\n")
             utils::flush.console()
          }
        }
      }
    }
  }
  return(invisible(out2))
}
