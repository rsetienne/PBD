#' @export

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
  outCRO = pbd_ML(brts = brts,initparsopt = initparsoptCR,idparsopt = 1:2,idparsfix = 3,parsfix = 10^8,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
  cat("Estimating parameters under the protracted birth-death model ...\n")
  outPBDO = pbd_ML(brts = brts,initparsopt = initparsoptPBD,idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
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
    outCR = pbd_ML(brtsCR[[mc]],initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = infinity,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
    outPBD1 = pbd_ML(brtsCR[[mc]],initparsopt = parsPBD[1:3],idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
    outPBD2 = pbd_ML(brtsCR[[mc]],initparsopt = c(parsCR + 0.05,length(brts) + 1000),idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
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
    outCR = pbd_ML(brtsPBD[[mc]],initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = infinity,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
    outPBD1 = pbd_ML(brtsPBD[[mc]],initparsopt = parsPBD[1:3],idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
    outPBD2 = pbd_ML(brtsPBD[[mc]],initparsopt = c(parsCR + 0.05,length(brts) + 1000),idparsopt = 1:3,idparsfix = NULL,parsfix = NULL,exteq = 1,parsfunc = parsfunc,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,methode = methode,n_low = n_low,n_up = n_up,tol = tol,maxiter = maxiter, optimmethod = optimmethod, verbose = verbose)
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
