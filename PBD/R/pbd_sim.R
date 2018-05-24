#' Function to simulate the protracted speciation process
#'
#' Simulating the protracted speciation process using the Doob-Gillespie
#' algorithm. This function differs from pbd_sim_cpp that 1) it does not
#' require that the speciation-initiation rate is the same for good and
#' incipient species, and 2) that it simulates the exact protracted speciation
#' process, and not the approximation made by the coalescent point process.
#' This function provides also the conversion to the approximation as output.
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b_1,
#' the speciation-initiation rate of good species \cr \code{pars[2]}
#' corresponds to la_1, the speciation-completion rate \cr \code{pars[3]}
#' corresponds to b_2, the speciation-initiation rate of incipient species \cr
#' \code{pars[4]} corresponds to mu_1, the extinction rate of good species \cr
#' \code{pars[5]} corresponds to mu_2, the extinction rate of incipient species
#' \cr
#' @param age Sets the age for the simulation
#' @param soc Sets whether this age is the stem (1) or crown (2) age
#' @param plotit Sets whether the various trees produced by the function should
#' be plotted or not
#' @param limitsize Sets a maximum to the number of incipient + good species
#' that are created during the simulation; if exceeded, the simulation is
#' aborted and removed.
#' @return \item{out}{ A list with the following elements: \cr \cr \code{tree}
#' is the tree of extant species in phylo format \cr \code{stree_random} is a
#' tree with one random sample per species in phylo format \cr
#' \code{stree_oldest} is a tree with the oldest sample per species in phylo
#' format \cr \code{stree_youngest} is a tree with the youngest sample per
#' species in phylo format \cr \code{L} is a matrix of all events in the
#' simulation where \cr - the first column is the incipient-level label of a
#' species \cr - the second column is the incipient-level label of the parent
#' of the species \cr - the third column is the time at which a species is born
#' as incipient species\cr - the fourth column is the time of
#' speciation-completion of the species \cr If the fourth element equals -1,
#' then the species is still incipient.  - the fifth column is the time of
#' extinction of the species \cr If the fifth element equals -1, then the
#' species is still extant.  - The sixth column is the species-level label of
#' the species \cr \code{sL_random} is a matrix like L but for
#' \code{stree_random} \cr \code{sL_oldest} is a matrix like L but for
#' \code{stree_oldest} \cr \code{sL_youngest} is a matrix like L but for
#' \code{stree_youngest} \cr \code{igtree.extinct} is the tree in simmap format
#' with incipient and good flags and including extinct species \cr
#' \code{igtree.extant} is the tree in simmap format with incipient and good
#' flags without extinct species \cr \code{recontree} is the reconstructed tree
#' in phylo format, reconstructed using the approximation in Lambert et al.
#' 2014 \cr \code{reconL} is the matrix corresponding to \code{recontree} \cr
#' \code{L0} is a matrix where the crown age is at 0; for internal use only \cr
#' }
#' @author Rampal S. Etienne
#' @seealso \code{\link{pbd_sim_cpp}}
#' @keywords models
#' @examples
#'  pbd_sim(c(0.2,1,0.2,0.1,0.1),15)
#' @export pbd_sim
pbd_sim = function(pars,age,soc = 2,plotit = FALSE, limitsize = Inf)
{
la1 = pars[1]
la2 = pars[2]
la3 = pars[3]
mu1 = pars[4]
mu2 = pars[5]

i = 1
while(i <= soc)
{
   t = 0
   if(i == 1)
   {
      id1 = 0
      id = id1 + 1
      Sid1 = 0
      Sid = 1
      sg = id
      si = NULL
      L = t(as.matrix(c(id,0,-1e-10,t,-1,1)))
   }
   if(i == 2)
   {
      id = id1 + 1
      Sid = Sid1
      sg = NULL
      si = -id
      L = t(as.matrix(c(id,1,t,-1,-1,1)))
   }

   Ng = length(sg)
   Ni = length(si)
   probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)
   denom = sum(probs)
   probs = probs/denom
   t = t - log(stats::runif(1))/denom

   while(t <= age)
   {
      event = DDD::sample2(1:5,size = 1,prob = probs)
      if(event == 1)
      {
          parent = as.numeric(DDD::sample2(sg,1))
          id = id + 1
          L = rbind(L,c(id,parent,t,-1,-1,L[abs(parent) - id1,6]))
          si = c(si,-id)
          Ni = Ni + 1
      }
      if(event == 2)
      {
          todie = as.numeric(DDD::sample2(sg,1))
          L[todie - id1,5] = t
          sg = sg[-which(sg == todie)]
          Ng = Ng - 1
      }
      if(event == 3)
      {
          tocomplete = abs(as.numeric(DDD::sample2(si,1)))
          L[tocomplete - id1,4] = t
          Sid = Sid + 1
          L[tocomplete - id1,6] = Sid
          sg = c(sg,tocomplete)
          si = si[-which(abs(si) == tocomplete)]
          Ng = Ng + 1
          Ni = Ni - 1
      }
      if(event == 4)
      {
          parent = as.numeric(DDD::sample2(si,1))
          id = id + 1
          L = rbind(L,c(id,parent,t,-1,-1,L[abs(parent) - id1,6]))
          si = c(si,-id)
          Ni = Ni + 1
      }
      if(event == 5)
      {
          todie = abs(as.numeric(DDD::sample2(si,1)))
          L[todie - id1,5] = t
          si = si[-which(abs(si) == todie)]
          Ni = Ni - 1
      }
      if(Ng + Ni > limitsize)
      {
        Ni = 0
        Ng = 0
      }
      probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)
      denom = sum(probs)
      probs = probs/denom
      t = t - log(stats::runif(1))/denom
   }
   if(i == 1)
   {
       if((Ng + Ni) > 0)
       {
           i = i + 1
           L1 = L
           id1 = id
           Sid1 = Sid
           si1 = si
           sg1 = sg
       }
   } else {
   if(i == 2)
   {
       if(checkgood(L,si,sg) == 1)
       {
           i = i + 1
           L2 = L
           si2 = si
           sg2 = sg
       }
   }}
}
L = L1
if(soc == 2)
{
    L = rbind(L1,L2)
}
L0 = L
absL = L
absL[,2] = abs(L[,2])
trans = NULL
igtree.extinct = phytools::read.simmap(text = detphy(absL,age,ig = T,dropextinct = F))
igtree.extant = phytools::read.simmap(text = detphy(absL,age,ig = T,dropextinct = T))
tree = ape::as.phylo(ape::read.tree(text = detphy(absL,age)))

# Random sampling
sL_random = sampletree(absL, age, samplemethod = "random")
stree_random = ape::as.phylo(ape::read.tree(text = detphy(sL_random, age)))
# Sampling the oldest
sL_oldest = sampletree(absL, age, samplemethod = "oldest")
stree_oldest = ape::as.phylo(ape::read.tree(text = detphy(sL_oldest, age)))
# Sampling the longest
sL_longest <- sampletree(absL, age, samplemethod = "longest")
stree_longest <- ape::as.phylo(ape::read.tree(text = detphy(sL_longest, age)))

sL_random[, 3:5][which(sL_random[, 3:5] == -1)] = age + 1
sL_random[, 3:5] = age - sL_random[, 3:5]
sL_random = sL_random[order(sL_random[, 1]), ]
sL_oldest[, 3:5][which(sL_oldest[, 3:5] == -1)] = age + 1
sL_oldest[, 3:5] = age - sL_oldest[, 3:5]
sL_oldest = sL_oldest[order(sL_oldest[, 1]), ]
sL_youngest[, 3:5][which(sL_youngest[, 3:5] == -1)] = age + 1
sL_youngest[, 3:5] = age - sL_youngest[, 3:5]
sL_youngest = sL_youngest[order(sL_youngest[, 1]), ]

#sL = sampletree(absL,age)
#stree = as.phylo(read.tree(text = detphy(sL,age)))
#sL[,3:5][which(sL[,3:5] == -1)] = age + 1
#sL[,3:5] = age - sL[,3:5]
#sL = sL[order(sL[,1]),]

#rbrts = sort(age - pbd_sim_step2a(L)[[1]])
reconL = pbd_reconstruct(L0)
recontree = ape::as.phylo(ape::read.tree(text = detphy(reconL,age)))
L[,3:5][which(L[,3:5] == -1)] = age + 1
L[,3:5] = age - L[,3:5]
L = L[order(L[,1]),]
#reducL = cbind(L[,3],abs(L[,2:1]),L[,5],L[,4],L[,6])
#fulltree = L2phylo2(reducL,dropextinct = F)
reconL[,3:5][which(reconL[,3:5] == -1)] = age + 1
reconL[,3:5] = age - reconL[,3:5]
reconL = reconL[order(reconL[,1]),]

if(plotit == TRUE)
{
   graphics::par(mfrow = c(3,3))
   cols = stats::setNames(c("gray","black"),c("i","g"))
   phytools::plotSimmap(igtree.extinct,colors = cols)
   phytools::plotSimmap(igtree.extant,colors = cols)
   graphics::plot(tree, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   graphics::plot(stree_random, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   graphics::plot(stree_oldest, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   graphics::plot(stree_youngest, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   graphics::plot(recontree, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   graphics::par(mfrow = c(1,1))
}

Ltreeslist = list(tree = tree, stree_random = stree_random, stree_oldest = stree_oldest, stree_youngest = stree_youngest, L = L, sL_random = sL_random, sL_oldest = sL_oldest, sL_youngest = sL_youngest, igtree.extinct = igtree.extinct, igtree.extant = igtree.extant, recontree = recontree, reconL = reconL, L0 = L0)
return(Ltreeslist)
}
