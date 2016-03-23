#' Likelihood maximization for protracted birth-death model of diversification
#'
#' @param brts branching times (positive, from present to past).
#' @note 'max(brts)' equals the crown age and
#'   'min(brts)' equals the most recent branching time
#' @param initparsopt the initial parameter values:
#'   initparsopt[1] = b (= la_1 in ER2012) = speciation initiation rate.
#'   initparsopt[2] = mu_1 (= mu_g in ER2012) = extinction rate of good species.
#'   initparsopt[3] = la_1 (= la_2 in ER2012) = speciation completion rate.
#'   initparsopt[4] = mu_2 (= mu_i in ER2012) = extinction rate of incipient species.
#' @param exteq incipient species have the same (1) or different (0) extinction rate as good species
#' @param parsfunc functions of parameters
#' @param res resolution of the method; res should be larger than the total number of species
#' @param missnumspec number of missing species
#' @param cond conditioning:
#'   cond = 0 conditioning on stem or clade age.
#'   cond = 1 conditioning on age and non-extinction of the phylogeny.
#'   cond = 2 conditioning on age and on number of extant taxa
#' @param btorph likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
#' @param soc stem (1) or crown (2) age
#' @param methode method of the numerical integration; see package deSolve for details
#' @param n_low lower bound on number of species (cond = 2)
#' @param n_up upper bound on number of species (cond = 2)
#' @param tol tolerance in optimization:
#    reltolx = relative tolerance of parameter values in optimization.
#    reltolf = relative tolerance of function value in optimization.
#    abstolx = absolute tolerance of parameter values in optimization
#' @param maxiter the maximum number of iterations in the optimization
#' @param optimmethod 'subplex' (current default) or 'simplex' (default of previous versions)
#' @param verbose create output while performing algorithm
#' @seealso bdd_loglik
#' @return Logarithm of the maximum likelihood estimate
#' @author Rampal S. Etienne
#' @examples
#'   pbd_ML(1:10,initparsopt = c(0.2,0.01,0.3), exteq = 1)
#' @export
pbd_ML = function(
  brts,
  initparsopt = c(0.2,0.1,1),
  idparsopt = 1:length(initparsopt),
  idparsfix = NULL,
  parsfix = NULL,
  exteq = 1,
  parsfunc = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), missnumspec = 0, cond = 1, btorph = 1, soc = 2, methode = "lsoda", n_low = 0, n_up = 0, tol = c(1E-6, 1E-6, 1E-6), maxiter = 1000 * round((1.25)^length(idparsopt)),
  optimmethod = 'subplex',
  verbose = TRUE)
{

  options(warn=-1)
  brts = sort(abs(as.numeric(brts)),decreasing = TRUE)
  if(is.numeric(brts) == FALSE)
  {
    cat("The branching times should be numeric.\n")
    out2 = data.frame(b = -1, mu_1 = -1, lambda_1 = -1,mu_2 = -1, loglik = -1, df = -1, conv = -1)
  } else {
    if(exteq == 1){ idexteq = 4 } else { idexteq = NULL }
    idpars = sort(c(idparsopt,idparsfix,idexteq))
    if((sum(idpars == (1:4)) != 4) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
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
      flush.console()
      initloglik = pbd_loglik_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,exteq = exteq,parsfunc = parsfunc,pars2 = pars2,brts = brts,missnumspec = missnumspec)
      if (verbose) { cat("The likelihood for the initial parameter values is",initloglik,"\n") }
      flush.console()
      if(initloglik == -Inf)
      {
        cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
        out2 = data.frame(b = -1, mu_1 = -1, lambda_1 = -1,mu_2 = -1, loglik = -1, df = -1, conv = -1)
      } else {
        if (verbose) { cat("Optimizing the likelihood - this may take a while.","\n") }
        flush.console()
        optimpars = c(tol,maxiter)
        #out = pbd_simplex(trparsopt,idparsopt,trparsfix,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
        out = optimizer(optimmethod = optimmethod,optimpars = optimpars,fun = pbd_loglik_choosepar,trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,exteq = exteq, parsfunc = parsfunc, pars2 = pars2,brts = brts, missnumspec = missnumspec)
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
          s1 = sprintf('Maximum likelihood parameter estimates: b: %f, mu_1: %f, lambda_1: %f, mu_2: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4])
          s2 = sprintf('Maximum loglikelihood: %f',ML)
          s3 = sprintf('The expected duration of speciation for these parameters is: %f',pbd_durspec_mean(c(MLpars1[1],MLpars1[3],MLpars1[4])))
          s4 = sprintf('The median duration of speciation for these parameters is: %f',pbd_durspec_quantile(c(MLpars1[1],MLpars1[3],MLpars1[4]),0.5))
          if (verbose) {
            cat("\n",s1,"\n",s2,"\n",s3,"\n",s4,"\n")
            flush.console()
          }
        }
      }
    }
  }
  invisible(out2)
}
