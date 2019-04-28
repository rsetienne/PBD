#' Bootstrap likelihood ratio test of protracted birth-death model of
#' diversification
#'
#' This function computes the maximum likelihood and the associated estimates
#' of the parameters of a protracted birth-death model of diversification for a
#' given set of phylogenetic branching times. It then performs a bootstrap
#' likelihood ratio test of the protracted birth-death (PBD) model against the
#' constant-rates (CR) birth-death model. Finally, it computes the power of
#' this test.
#'
#' The output is a list with 3 elements:
#'
#' @param brts A set of branching times of a phylogeny, all positive
#' @param initparsoptPBD The initial values of the parameters that must be
#' optimized for the protracted birth-death (PBD) model: b, mu and lambda
#' @param initparsoptCR The initial values of the parameters that must be
#' optimized for the constant-rates (CR) model: b and mu
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @param outputfilename The name (and location) of the file where the output
#' will be saved. Default is no save.
#' @param seed The seed for the pseudo random number generator for simulating
#' the bootstrap data
#' @param endmc The number of bootstraps
#' @param alpha The significance level of the test
#' @param plotit Boolean to plot results or not
#' @param parsfunc Specifies functions how the rates depend on time, default
#' functions are constant functions
#' @param cond Conditioning: \cr
#' cond == 0 : conditioning on stem or crown age\cr
#' cond == 1 : conditioning on stem or crown age and non-extinction of the
#' phylogeny \cr
#' cond == 2 : conditioning on stem or crown age and on the total
#' number of extant taxa (including missing species) \cr
#' cond == 3 : conditioning on the total number of extant taxa (including
#' missing species)\cr
#' Note: cond == 3 assumes a uniform prior on stem age, as is the standard
#' in constant-rate birth-death models, see e.g. D. Aldous & L. Popovic 2004.
#' Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J. Theor. Biol. 261:
#' 58-66.
#' @param btorph Sets whether the likelihood is for the branching times (0) or
#' the phylogeny (1)
#' @param soc Sets whether stem or crown age should be used (1 or 2)
#' @param methode The numerical method used to solve the master equation, such
#' as 'lsoda' or 'ode45'.
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
#' default is 'subplex'. Alternative is 'simplex'.
#' @param num_cycles Number of cycles of the optimization (default is 1).
#' @param verbose if TRUE, explanatory text will be shown
#' @return \item{brtsCR}{a list of sets of branching times generated under the
#' constant-rates model using the ML parameters under the CR model}
#' \item{brtsDD}{a list of sets of branching times generated under the
#' protracted birth-death model using the ML parameters under the PBD model}
#' \item{out}{a dataframe with the parameter estimates and maximum likelihoods
#' for protracted birth-death and constant-rates models
#' \code{$model} - the model used to generate the data. 0 = unknown (for real
#' data), 1 = CR, 2 = PBD \cr
#' \code{$mc} - the simulation number for each model \cr
#' \code{$b_CR} - peciation rate estimated under CR \cr
#' \code{$mu_CR} - extinction rate estimated under CR \cr
#' \code{$LL_CR} - maximum likelihood estimated under CR\cr
#' \code{$conv_CR} - convergence code for likelihood optimization; conv = 0
#' means convergence \cr
#' \code{$b_PBD1} - speciation-initation rate estimated
#' under PBD for first set of initial values\cr
#' \code{$mu_PB1} - extinction
#' rate estimated under DD for first set of initial values \cr
#' \code{$lambda_PB1} - speciation-completion rate estimated under PBD for
#' first set of initial values \cr
#' \code{$LL_PBD1} - maximum likelihood estimated under DD for first set of
#' initial values \cr \code{$conv_PBD1} -
#' convergence code for likelihood optimization for first set of initial
#' values; conv = 0 means convergence \cr
#' \code{$b_PBD2} - speciation-initation
#' rate estimated under PBD for second set of initial values\cr
#' \code{$mu_PB2} - extinction rate estimated under DD for second set of
#' initial values \cr
#' \code{$lambda_PB2} - speciation-completion rate estimated under PBD for
#' second set of initial values \cr
#' \code{$LL_PBD2} - maximum likelihood
#' estimated under DD for second set of initial values \cr
#' \code{$conv_PBD2} - convergence code for likelihood optimization for
#' econd set of initial values; conv = 0 means convergence \cr
#' \code{$LR} - likelihood ratio between DD and CR }
#' \item{pvalue}{p-value of the test}
#' \item{LRalpha}{Likelihood
#' ratio at the signifiance level alpha}
#' \item{poweroftest}{power of the test for significance level alpha}
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_loglik}}, \code{\link{pbd_ML}}
#' @references - Etienne, R.S. et al. 2016. Meth. Ecol. Evol. 7: 1092-1099,
#' doi: 10.1111/2041-210X.12565 \cr
#' @keywords models
#' @export pbd_LR
pbd_LR = function(
  brts,
  initparsoptPBD,
  initparsoptCR,
  missnumspec,
  outputfilename = NULL,
  seed = 42,
  endmc = 1000,
  alpha = 0.05,
  plotit = TRUE,
  parsfunc = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}),
  cond = 1,
  btorph = 1,
  soc = 2,
  methode = 'lsoda',
  n_low = 0,
  n_up = 0,
  tol = c(1E-6,1E-6,1E-6),
  maxiter = 2000,
  optimmethod = 'subplex',
  num_cycles = num_cycles,
  verbose = FALSE
)
{
  if(!is.null(seed))
  {
    set.seed(DDD::roundn(seed))
  }
  infinity = 10^10
  age = max(brts)
  cat("Estimating parameters under the constant-rate model ...\n")
  outCRO = pbd_ML(brts = brts,initparsopt = initparsoptCR,idparsopt = 1:2,idparsfix = 3,parsfix = 10^8,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
  cat("Estimating parameters under the protracted birth-death model ...\n")
  outPBDO = pbd_ML(brts = brts,initparsopt = initparsoptPBD,idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
  LRO = outPBDO$loglik - outCRO$loglik
  out = cbind(NA,NA,outCRO,outPBDO,NA,NA,NA,NA,NA,LRO)
  out = out[,-c(5,6,8,13,15)]
  newnames = c("model","mc","b_CR","mu_CR","LL_CR","conv_CR","b_PBD1","mu_PBD1","lambda_PBD1","LL_PBD1","conv_PBD1","b_PBD2","mu_PBD2","lambda_PBD2","LL_PBD2","conv_PBD2","LR")
  names(out) = newnames
  if(!is.null(outputfilename))
  {
    save(seed,brts,out,file = outputfilename)
  }
  parsCR = as.numeric(outCRO[1:2])
  parsPBD = as.numeric(outPBDO[1:4])
  brtsCR = list()
  brtsPBD = list()
  cat('Simulating trees under CR and PBD models ...\n')
  for(mc in 1:endmc)
  {
    cat('Simulation pair',mc,'\n')
    brtsCR[[mc]] = pbd_sim_cpp(pars = c(parsCR,infinity,parsCR[2]),age = age,parsf = parsfunc,soc = soc,plotltt = 0,methode = methode)
    brtsPBD[[mc]] = pbd_sim_cpp(pars = parsPBD,age = age,parsf = parsfunc,soc = soc,plotltt = 0,methode = methode)
    utils::flush.console()
  }
  if(!is.null(outputfilename))
  {
    save(seed,brts,out,brtsCR,brtsPBD,file = outputfilename)
  }
  cat('Performing bootstrap to determine critical LR ...\n')
  for(mc in 1:endmc)
  {
    cat('Analyzing simulation:',mc,'\n')
    utils::flush.console()
    outCR = pbd_ML(brtsCR[[mc]],initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = infinity,exteq = 1,parsfunc = parsfunc,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
    outPBD1 = pbd_ML(brtsCR[[mc]],initparsopt = parsPBD[1:3],idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
    outPBD2 = pbd_ML(brtsCR[[mc]],initparsopt = c(parsCR + 0.05,length(brts) + 1000),idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
    if(outPBD1$conv == -1 & outPBD2$conv == -1)
    {
      maxLLPBD = outCR$loglik
    } else if(outPBD1$conv == -1 & outPBD2$conv != -1)
    {
      maxLLPBD = outPBD2$loglik
    } else if(outPBD1$conv != -1 & outPBD2$conv == -1)
    {
      maxLLPBD = outPBD1$loglik
    } else {
      maxLLPBD = max(outPBD1$loglik,outPBD2$loglik)
    }
    LR = pmax(0,maxLLPBD - outCR$loglik)
    outff = cbind(1,mc,outCR,outPBD1,outPBD2,LR)
    outff = outff[,-c(5,6,8,13,15,20,22)]
    names(outff) = newnames
    out = rbind(out,outff)
    if(!is.null(outputfilename))
    {
      save(seed,brts,out,brtsCR,brtsPBD,file = outputfilename)
    }
  }
  opt = rep(0,endmc)
  cat('Performing bootstrap to determine power ...\n')
  for(mc in 1:endmc)
  {
    cat('Analyzing simulation:',mc,'\n')
    utils::flush.console()
    outCR = pbd_ML(brtsPBD[[mc]],initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = infinity,exteq = 1,parsfunc = parsfunc,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod,num_cycles = num_cycles, verbose = verbose)
    outPBD1 = pbd_ML(brtsPBD[[mc]],initparsopt = parsPBD[1:3],idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
    outPBD2 = pbd_ML(brtsPBD[[mc]],initparsopt = c(parsCR + 0.05,length(brts) + 1000),idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, num_cycles = num_cycles, verbose = verbose)
    if(outPBD1$conv == -1 & outPBD2$conv == -1)
    {
      maxLLPBD = outCR$loglik
      opt[mc] = 1
    } else if(outPBD1$conv != -1 & outPBD2$conv == -1)
    {
      maxLLPBD = outPBD1$loglik
      opt[mc] = 2
    } else if(outPBD1$conv == -1 & outPBD2$conv != -1)
    {
      maxLLPBD = outPBD2$loglik
      opt[mc] = 3
    } else {
      maxLLPBD = max(outPBD1$loglik,outPBD2$loglik)
      opt[mc] = 1 + min(which(c(outPBD1$loglik,outPBD2$loglik) == maxLLPBD))
    }
    LR = pmax(0,maxLLPBD - outCR$loglik)
    outff = cbind(1,mc,outCR,outPBD1,outPBD2,LR)
    outff = outff[,-c(5,6,8,13,15,20,22)]
    names(outff) = newnames
    out = rbind(out,outff)
    if(!is.null(outputfilename))
    {
      save(seed,brts,out,opt,brtsCR,brtsPBD,file = outputfilename)
    }
  }
  inverse_quantile = function(samples,x)
  {
    samplessort = sort(samples)
    pup = which(samplessort > x)
    if(length(pup) > 0)
    {
      if(length(pup) < length(samplessort))
      {
        pup = min(pup)
        invquant = (pup + (x - samplessort[pup])/(samplessort[pup - 1] - samplessort[pup]))/length(samples)
      } else {
        invquant = 0
      }
    } else {
      invquant = 1
    }
    return(invquant)
  }
  funpvalue = function(samples,x)
  {
    samplessort = sort(samples)
    pup = which(samplessort > x)
    pvalue = (length(pup) + 1)/ (length(samples) + 1)
    return(pvalue)
  }
  funpoweroftest = function(samples,x)
  {
    samplessort = sort(samples)
    pup = which(samplessort > x)
    poweroftest = length(pup)/(length(samples) + 1)
    return(poweroftest)
  }
  #pvalue = 1 - inverse_quantile(out$LR[2:(endmc + 1)],out$LR[1])
  pvalue = funpvalue(out$LR[2:(endmc + 1)],out$LR[1])
  LRalpha = as.numeric(stats::quantile(out$LR[2:(endmc + 1)],1 - alpha,type = 4))
  #poweroftest = 1 - inverse_quantile(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  poweroftest = funpoweroftest(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  if(plotit == TRUE)
  {
    try(grDevices::dev.off())
    try(grDevices::dev.off())
    pdffilename = paste(getwd(),'/LR.pdf',sep = '')
    grDevices::pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
    al = 0.03
    alw = 2
    alw2 = 1.7
    aa = 45
    graphics::par(mfrow = c(2,2),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
    graphics::hist(out$LR[2:(1 + endmc)],main = 'Distribution of LLR under CR',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(out$LR[1:(endmc + 1)])))
    graphics::arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
    graphics::arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')
    graphics::box()
    graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", col = "white", col.lab = 'white', col.axis = 'white')
    graphics::hist(out$LR[(endmc + 2):(1 + 2 * endmc)], main = 'Distribution of LLR under PBD',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30)
    graphics::box()
    graphics::arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
    graphics::arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')

    graphics::par(mfrow = c(2,3),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
    b = out$b_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$b_PBD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$b_PBD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
    mu = out$mu_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$mu_PBD1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$mu_PBD2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
    lambda = 1E+120 * (opt == 1) + pmin(1E+120,out$lambda_PBD1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$lambda_PBD2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
    graphics::hist(b,main = NULL, xlab = expression(b), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(b)))
    graphics::arrows(out$b_PBD1[1],-1E+120, x1 = out$b_PBD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
    graphics::box()
    graphics::hist(mu,main = NULL, xlab = expression(mu), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(mu)))
    graphics::arrows(out$mu_PBD1[1],-1E+120, x1 = out$mu_PBD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
    graphics::box()
    graphics::hist(lambda,main = NULL, xlab = expression(lambda), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(lambda)))
    graphics::arrows(out$lambda_PBD1[1],-1E+120, x1 = out$lambda_PBD1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
    graphics::box()
    try(grDevices::dev.off())
    try(grDevices::dev.off())
    os = .Platform$OS.type
    if(os == "windows")
    {
      shell.exec(pdffilename)
    }
    if(os == "unix")
    {
      system(paste("open",pdffilename,sep = " "))
    }
  }
  if(!is.null(outputfilename))
  {
    save(seed,brts,out,opt,brtsCR,brtsPBD,pvalue,LRalpha,poweroftest,file = outputfilename)
  }
  return(list(brtsCR = brtsCR,brtsPBD = brtsPBD,out = out,pvalue = pvalue,LRalpha = LRalpha,poweroftest = poweroftest))
}
