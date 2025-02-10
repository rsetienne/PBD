---
title: "pbd_ML demo"
author: "Richel Bilderbeek & Rampal S. Etienne"
date: "2017-04-28"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{pbd_ML demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This document gives a demonstration how to use 
the package to obtain a maximum-likelihood estimate
of the protracted birth-death speciation model.

First thing is to load the PBD package itself:



```r
rm(list = ls())
library(PBD)
```

We will also need ape for `branching.times`:


```r
library(ape)
```


Here we simulate a tree with known parameters:


```r
seed <- 43
set.seed(seed)
b_1 <- 0.3 # speciation-initiation rate of good species
la_1 <- 0.2 # speciation-completion rate
b_2 <- b_1 # the speciation-initiation rate of incipient species
mu_1 <- 0.1 #  extinction rate of good species
mu_2 <- mu_1 # extinction rate of incipient species 
pars <- c(b_1, la_1, b_2, mu_1, mu_2)
age <- 15 # the age for the simulation 
phylogenies <- pbd_sim(pars = pars, age = age)
plot(phylogenies$$recontree)
plot(phylogenies$igtree.extant)
plot(phylogenies$tree)
names(phylgenies
```

```
## Error: <text>:11:18: unexpected '$'
## 10: phylogenies <- pbd_sim(pars = pars, age = age)
## 11: plot(phylogenies$$
##                      ^
```

Now we try to recover the parameters by maximum likelihood estimation:



```r
brts <- branching.times(phylogeny)  # branching times
```

```
## Error in inherits(phy, "phylo"): object 'phylogeny' not found
```

```r
init_b <- 0.2  # speciation-initiation rate
init_mu_1 <- 0.05  # extinction rate of good species
init_la_1 <- 0.3 # speciation-completion rate
#init_mu_2 <- 0.05  # extinction rate of incipient species  # not used

# The initial values of the parameters that must be optimized
initparsopt <- c(init_b, init_mu_1, init_la_1)

# The extinction rates between incipient and good species are equal
exteq <- TRUE

# The first element of the branching times is the crown age (and not the stem age)
soc <- 2

# Conditioning on non-extinction of the phylogeny
# as I actively selected for a nice phylogeny
cond <- 1

# Give the likelihood of the phylogeny (instead of the likelihood of the branching times)
btorph <- 1
```

Maximum likelihood estimation can now be performed:


```r
r <- pbd_ML(
  brts = brts,
  initparsopt = initparsopt, 
  exteq = exteq,
  soc = soc, 
  cond = cond,
  btorph = btorph,
  verbose = FALSE
)
```

```
## Error in sort(abs(as.numeric(brts)), decreasing = TRUE): object 'brts' not found
```

The ML parameter estimates are:


```r
knitr::kable(r)
```

```
## Error in inherits(x, "list"): object 'r' not found
```

Comparing the known true value with the recovered values:


```r
loglik_true <- PBD::pbd_loglik(pars, brts = brts)
```

```
## Error in PBD::pbd_loglik(pars, brts = brts): object 'pars' not found
```

```r
df <- as.data.frame(r)
```

```
## Error in as.data.frame(r): object 'r' not found
```

```r
df <- rbind(df, c(b_1, mu_1, la_1, mu_2, loglik_true, NA, NA))
```

```
## Error in rbind(df, c(b_1, mu_1, la_1, mu_2, loglik_true, NA, NA)): object 'b_1' not found
```

```r
row.names(df) <- c("ML", "true")
```

```
## Error in `rownames<-`(x, value): attempt to set 'rownames' on an object with no dimensions
```

```r
knitr::kable(df)
```

```
## Error in as.data.frame.default(x): cannot coerce class ""function"" to a data.frame
```

Ideally, all parameter columns should have the same values.

To test for the uncertainty of our ML estimate, we can do a parametric bootstrap.

The function `pbd_bootstrap` does:

 * First do a ML estimate
 * Run a simulation with those estimates, then recover these estimates by ML estimation


```r
endmc <- 10 # Sets the number of simulations for the bootstrap

b <- pbd_bootstrap(
  brts = brts,
  initparsopt = initparsopt, 
  exteq = exteq,
  soc = soc, 
  cond = cond,
  btorph = btorph,
  plotltt = FALSE,
  endmc = endmc,
  seed = seed
)
```

```
## Error in set.seed(seed): object 'seed' not found
```

```r
knitr::kable(b[[3]])
```

```
## Error in knitr::kable(b[[3]]): object 'b' not found
```

From the bootstrap analysis, we get 

 * Again the ML estimate
 * The ML estimates for simulations run with those estimates

Putting this in a table:


```r
dg <- rbind(df, 
  list(
    b = b[[1]]$b, 
    mu_1 = b[[1]]$mu_1, 
    lambda_1 = b[[1]]$lambda_1, 
    mu_2 = b[[1]]$mu_2,
    loglik = b[[1]]$loglik,
    df = b[[1]]$df,
    conv = b[[1]]$conv
  ),
  list(
    b = b[[3]]$b, 
    mu_1 = b[[3]]$mu_1, 
    lambda_1 = b[[3]]$lambda_1, 
    mu_2 = b[[3]]$mu_2,
    loglik = b[[3]]$loglik,
    df = b[[3]]$df,
    conv = b[[3]]$conv
  )
)
```

```
## Error in rbind(df, list(b = b[[1]]$b, mu_1 = b[[1]]$mu_1, lambda_1 = b[[1]]$lambda_1, : object 'b' not found
```

```r
dg
```

```
## Error in eval(expr, envir, enclos): object 'dg' not found
```

```r
row.names(dg) <- c("ML", "true", "ML2", paste("BS", 1:endmc, sep = ""))
```

```
## Error in row.names(dg) <- c("ML", "true", "ML2", paste("BS", 1:endmc, : object 'dg' not found
```

```r
knitr::kable(dg)
```

```
## Error in knitr::kable(dg): object 'dg' not found
```

We expect rows ML and ML2 to be identical. Their values are
indeed very similar.

We can calculate the loglikelihood for 


```r
ml_b <- b[[1]]$b
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
ml_mu_1 <- b[[1]]$mu_1
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
ml_la_1 <- b[[1]]$lambda_1
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
ml_mu_2 <- b[[1]]$mu_2
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
ml_pars1 <- c(ml_b, ml_mu_1, ml_la_1, ml_mu_2)
```

```
## Error in eval(expr, envir, enclos): object 'ml_b' not found
```

```r
ml_pars2 <- c(cond, btorph, soc, 0, "lsoda")

l <- pbd_loglik(
  pars1 = ml_pars1,
  pars2 = ml_pars2,
  brts = brts
)
```

```
## Error in pbd_loglik(pars1 = ml_pars1, pars2 = ml_pars2, brts = brts): object 'ml_pars1' not found
```

```r
print(l)
```

```
## Error in print(l): object 'l' not found
```


```
# Create .md, .html, and .pdf files
setwd(paste(getwd(), "vignettes", sep = "/"))
knit("PBD_ML_demo.Rmd")
markdownToHTML('PBD_ML_demo.md', 'PBD_ML_demo.html', options=c("use_xhml"))
system("pandoc -s PBD_ML_demo.html -o PBD_ML_demo.pdf")
```
