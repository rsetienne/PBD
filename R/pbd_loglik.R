logcondfun = function(nn,socsoc,yy)
{
  if(socsoc == 1)
  {
    nfac = 0
  } else
  {
    nfac = log(nn - 1)
  }
  pnc = nfac + (nn - socsoc) * log(1 - yy)
  return(pnc)
}

#' Loglikelihood for protracted birth-death model of diversification
#'
#' pbd_loglik computes the loglikelihood of the parameters of the protracted
#' speciation model given a set of branching times and number of missing
#' species
#'
#'
#' @param pars1 Vector of parameters: \cr \cr \code{pars1[1]} corresponds to b
#' (= la_1 in Etienne & Rosindell R2012) = speciation initiation rate \cr
#' \code{pars1[2]} corresponds to mu_1 (= mu_g in Etienne & Rosindell 2012) =
#' extinction rate of good species \cr \code{pars1[3]} corresponds to la_1 (=
#' la_2 in Etienne & Rosindell 2012) = speciation completion rate \cr
#' \code{pars1[4]} corresponds to mu_2 (= mu_i in ER2012) = extinction rate of
#' incipient species \cr When rates depend on time this time dependence should
#' be specified in pars1f and pars1 then becomes the parameters used in pars1f
#' \cr \cr
#' @param pars1f Vector of functions how the rates depend on time, default
#' functions are constant functions of the parameters in pars1: \cr \cr
#' \code{pars1f[1]} corresponds to time-dependence of b (= la_1 in Etienne &
#' Rosindell R2012) = speciation initiation rate \cr \code{pars1f[2]}
#' corresponds to time-dependence of mu_1 (= mu_g in Etienne & Rosindell 2012)
#' = extinction rate of good species \cr \code{pars1f[3]} corresponds to
#' tiem-dependence of la_1 (= la_2 in Etienne & Rosindell 2012) = speciation
#' completion rate \cr \code{pars1f[4]} corresponds to time-dependence of mu_2
#' (= mu_i in ER2012) = extinction rate of incipient species \cr \cr
#' @param pars2 Vector of model settings: \cr \cr \code{pars2[1]} set the
#' conditioning on non-extinction of the clade (1) or not (0) \cr \cr
#' \code{pars2[2]} sets whether the likelihood is for the branching times (0)
#' or the phylogeny (1) \cr \cr \code{pars2[3]} sets whether the first element
#' of the branching times is the stem (1) or the crown (2) age \cr \cr
#' \code{pars2[4]} sets whether the parameters and likelihood should be shown
#' on screen (1) or not (0) \cr \cr \code{pars2[5]} sets which method should be
#' used in the ode-solver. Default is 'lsoda'. See package deSolve for details.
#' \cr \cr \code{pars2[6]}Sets the lower bound of the number of species on
#' which conditioning should be done when cond = 2. Set this to 0 when
#' conditioning should be done on precisely the number of species (default)\cr
#' \cr \code{pars2[7]}Sets the upper bound of the number of species on which
#' conditioning should be done when cond = 2. Set this to 0 when conditioning
#' should be done on precisely the number of species (default)\cr \cr
#' @param brts A set of branching times of a phylogeny, all positive
#' @param missnumspec The number of species that are in the clade but missing
#' in the phylogeny
#' @return The loglikelihood
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_ML}}
#' @keywords models
#' @examples
#'  pbd_loglik(pars1 = c(0.2,0.1,1,0.1), pars2 = c(1,1,2,0,"lsoda"),brts = 1:10)
#' @export pbd_loglik
pbd_loglik = function(
  pars1,pars1f = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}),
  pars2 = c(1,1,2,1,"lsoda",0,0),
  brts,
  missnumspec = 0
)
{
# pbd_loglik computes the loglikelihood of the protracted speciation model given a set of branching times and data

# pars1 contains model parameters
# In the simplest case where rates do not depend on time, we have
# - pars1[1] = b (= la_1 in ER2012) = speciation initiation rate
# - pars1[2] = mu_1 (= mu_g in ER2012) = extinction rate of good species
# - pars1[3] = la_1 (= la_2 in ER2012) = speciation completion rate
# - pars1[4] = mu_2 (= mu_i in ER2012) = extinction rate of incipient species
# When rates depend on time this time dependence should be specified in pars1f and pars1 then become the parameters used in pars1f
# pars1f contains the functions how the rates depend on time, default functions are constant functions of the parameters in pars1
# pars2 contains settings
# - pars2[1] = cond = conditioning on age (0), age and non-extinction of the clade (1) or age and number of extant taxa (2)
# - pars2[2] = btorph = likelihood for branching times (0) or phylogeny (1)
# - pars2[3] = soc = stem (1) or crown (2) age
# - pars2[4] = printing of parameters and likelihood (1) or not (0)
# - pars2[5] = method of the numerical integration; see package deSolve for details
# brts = set of branching times
# missnumspec = the number of species that belong to the same clade but are not in the phylogeny

# Example: pbd_loglik(pars1 = c(0.1,0.05,1,0.05), brts = 1:10, missnumspec = 4)

  pars1 = c(pars1f,pars1)

  brts = sort(abs(brts))
  abstol = 1e-16
  reltol = 1e-10
  b = pars1[[1]](brts,as.numeric(pars1[5:length(pars1)]))
  methode = pars2[5]
  cond = as.numeric(pars2[1])
  btorph = as.numeric(pars2[2])
  soc = as.numeric(pars2[3])
  S = length(brts) + (soc - 1)
  m = missnumspec

  probs = c(1,1,0,0)
  y = deSolve::ode(probs,c(0,brts),pbd_loglik_rhs,c(pars1),rtol = reltol,atol = abstol,method = methode)
  if(dim(y)[1] < length(brts) + 1) { return(-Inf) }

  loglik = (btorph == 0) * lgamma(S) +
    (cond == 0) * soc * (log(y[length(brts) + 1,2]) + log(1 - y[length(brts) + 1,3])) +
    (cond == 1) * soc * log(y[length(brts) + 1,2])
  if(length(brts) > 1)
  {
    loglik = loglik + sum(log(b) + log(y[2:length(brts),2]) + log(1 - y[2:length(brts),3]))
  }
  if(cond == 2)
  {
    n_l = as.numeric(pars2[6])
    n_u = as.numeric(pars2[7])
    if(n_l == 0 & n_u == 0)
    {
      n_l = S + m
      n_u = S + m
    } else if(n_l > (S + m) | n_u < (S + m) | n_u < n_l)
    {
      cat('Lower or upper boundary not possible.\n')
      return(-Inf)
    }
    if(n_l == soc & n_u == Inf)
    {
      logcond = -soc * log(y[length(brts) + 1,2])
    } else if(n_u == Inf)
    {
      n_u = n_l - 1
      n_l = soc
      logcond = log(y[length(brts) + 1,2]^(-2) - sum(exp(logcondfun(n_l:n_u,soc,y[(length(brts) + 1),2]))))
      if(logcond == -Inf)
      {
        cat('Catastrophic cancellation encountered. Trying now with 10000 as upper bound.\n')
        logcond = log(sum(exp(logcondfun((n_u + 1):10000,soc,y[(length(brts) + 1),2]))))
      }
    } else
    {
      logcond = log(sum(exp(logcondfun(n_l:n_u,soc,y[(length(brts) + 1),2]))))
    }
    loglik = loglik - logcond
  }
  if(m > 0)
  {
    if(soc == 1)
    {
      y2 = as.numeric(c(1 - y[2:(length(brts) + 1),2]))
    }
    if(soc == 2)
    {
      y2 = as.numeric(c(1 - y[2:(length(brts) + 1),2],1 - y[length(brts) + 1,2]))
    }
    x = rep(0,m + 1)
    x[1] = 1
    for(j in 1:S)
    {
      #x = convolve(x,rev((1:(m + 1)) * (y2[j]^(0:m))),type = 'open')[1:(m + 1)]
      x = DDD::conv(x,(1:(m + 1)) * (y2[j]^(0:m)))[1:(m+1)]
    }
    loglik = loglik + lgamma(S + 1) + lgamma(m + 1) - lgamma(S + m + 1) + log(x[m + 1])
  }

  if(as.numeric(pars2[4]) == 1)
  {
    pastetxt = paste('Parameters:',pars1[[5]][1],sep = ' ')
    for(i in 6:length(pars1))
    {
      pastetxt = paste(pastetxt,pars1[[i]][1],sep = ', ')
    }
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(pastetxt,s2,"\n",sep = "")
    utils::flush.console()
  }

  return(as.numeric(loglik))
}
