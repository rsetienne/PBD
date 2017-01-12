resample<- function (resample_L, age, plotit = FALSE) 
{
  
  L=resample_L # resample_L is the tree to resample. Must be the "L" element of the list returned by pbd_sim function
  L0 = L
  absL = L
  absL[, 2] = abs(L[, 2])
  trans = NULL
  
  sL_random = sampletree(absL, age, samplemethod = "random")
  stree_random = as.phylo(read.tree(text = detphy(sL_random, 
                                                  age)))
  sL_oldest = sampletree(absL, age, samplemethod = "oldest")
  stree_oldest = as.phylo(read.tree(text = detphy(sL_oldest, 
                                                  age)))
  sL_youngest = sampletree(absL, age, samplemethod = "youngest")
  stree_youngest = as.phylo(read.tree(text = detphy(sL_youngest, 
                                                    age)))
  sL_random[, 3:5][which(sL_random[, 3:5] == -1)] = age + 1
  sL_random[, 3:5] = age - sL_random[, 3:5]
  sL_random = sL_random[order(sL_random[, 1]), ]
  
  sL_oldest[, 3:5][which(sL_oldest[, 3:5] == -1)] = age + 1
  sL_oldest[, 3:5] = age - sL_oldest[, 3:5]
  sL_oldest = sL_oldest[order(sL_oldest[, 1]), ]
  
  sL_youngest[, 3:5][which(sL_youngest[, 3:5] == -1)] = age + 
    1
  sL_youngest[, 3:5] = age - sL_youngest[, 3:5]
  sL_youngest = sL_youngest[order(sL_youngest[, 1]), ]
 
 
  if (plotit == TRUE) {
    par(mfrow = c(1, 3))
    cols = setNames(c("gray", "black"), c("i", "g"))
   
    plot(stree_random, edge.width = 2, font = 1, label.offset = 0.1, 
         cex = 1)
    plot(stree_oldest, edge.width = 2, font = 1, label.offset = 0.1, 
         cex = 1)
    plot(stree_youngest, edge.width = 2, font = 1, label.offset = 0.1, 
         cex = 1)
   
    par(mfrow = c(1, 1))
  }
  Ltreeslist = list(stree_random = stree_random, 
                    stree_oldest = stree_oldest, stree_youngest = stree_youngest, 
                    L = L, sL_random = sL_random, sL_oldest = sL_oldest, 
                    sL_youngest = sL_youngest, L0 = L0)
  return(Ltreeslist)
}
