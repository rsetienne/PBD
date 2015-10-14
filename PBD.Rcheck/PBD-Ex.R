pkgname <- "PBD"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "PBD-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('PBD')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PBD-package")
### * PBD-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PBD-package
### Title: Protracted birth-death model of diversification
### Aliases: PBD-package PBD
### Keywords: package

### ** Examples

pbd_ML(1:10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PBD-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_ML")
### * pbd_ML

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_ML
### Title: Maximization of loglikelihood under protracted birth-death model
###   of diversification
### Aliases: pbd_ML
### Keywords: models

### ** Examples

pbd_ML(1:10,initparsopt = c(0.2,0.01,0.3), exteq = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_ML", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_bootstrap")
### * pbd_bootstrap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_bootstrap
### Title: Bootstrap analysis under protracted birth-death model of
###   diversification
### Aliases: pbd_bootstrap
### Keywords: models

### ** Examples

pbd_bootstrap(1:10,endmc = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_bootstrap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_brts_density")
### * pbd_brts_density

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_brts_density
### Title: Node depth probbaility density for protracted birth-death model
###   of diversification
### Aliases: pbd_brts_density
### Keywords: models

### ** Examples
 pbd_brts_density(pars1 = c(0.2,0.1,1,0.1), methode = "lsoda",brts = 1:10) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_brts_density", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_cumdensity")
### * pbd_durspec_cumdensity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_cumdensity
### Title: Cumulative density of duration of speciation under protracted
###   birth-death model of diversification
### Aliases: pbd_durspec_cumdensity
### Keywords: models

### ** Examples
 pbd_durspec_cumdensity(pars = c(0.5,0.3,0.1),3) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_cumdensity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_density")
### * pbd_durspec_density

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_density
### Title: Probability density for duration of speciation under protracted
###   birth-death model of diversification
### Aliases: pbd_durspec_density
### Keywords: models

### ** Examples
 pbd_durspec_density(pars = c(0.5,0.3,0.1), tau = 1) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_density", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_mean")
### * pbd_durspec_mean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_mean
### Title: Mean duration of speciation under protracted birth-death model
###   of diversification
### Aliases: pbd_durspec_mean
### Keywords: models

### ** Examples
 pbd_durspec_mean(pars = c(0.5,0.3,0.1)) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_mean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_mode")
### * pbd_durspec_mode

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_mode
### Title: mode of the duration of speciation under protracted birth-death
###   model of diversification
### Aliases: pbd_durspec_mode
### Keywords: models

### ** Examples
 pbd_durspec_mode(pars = c(0.5,0.3,0.1)) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_mode", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_moment")
### * pbd_durspec_moment

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_moment
### Title: Moments of duration of speciation under protracted birth-death
###   model of diversification
### Aliases: pbd_durspec_moment
### Keywords: models

### ** Examples
 pbd_durspec_moment(pars = c(0.5,0.3,0.1),2) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_moment", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_quantile")
### * pbd_durspec_quantile

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_quantile
### Title: Quantiles of duration of speciation under protracted birth-death
###   model of diversification
### Aliases: pbd_durspec_quantile
### Keywords: models

### ** Examples
 pbd_durspec_quantile(pars = c(0.5,0.3,0.1),0.5) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_quantile", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_durspec_var")
### * pbd_durspec_var

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_durspec_var
### Title: Variance in duration of speciation under protracted birth-death
###   model of diversification
### Aliases: pbd_durspec_var
### Keywords: models

### ** Examples
 pbd_durspec_var(pars = c(0.5,0.3,0.1)) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_durspec_var", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_loglik")
### * pbd_loglik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_loglik
### Title: Loglikelihood for protracted birth-death model of
###   diversification
### Aliases: pbd_loglik
### Keywords: models

### ** Examples
 pbd_loglik(pars1 = c(0.2,0.1,1,0.1), pars2 = c(1,1,2,0,"lsoda"),brts = 1:10) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_loglik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_sim")
### * pbd_sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_sim
### Title: Function to simulate the protracted speciation process
### Aliases: pbd_sim
### Keywords: models

### ** Examples
 pbd_sim(c(0.2,1,0.2,0.1,0.1),15) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pbd_sim_cpp")
### * pbd_sim_cpp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pbd_sim_cpp
### Title: Function to simulate the approximate protracted speciation
###   process
### Aliases: pbd_sim_cpp
### Keywords: models

### ** Examples
 pbd_sim_cpp(pars = c(0.2,1,0.2,0.1),age = 15) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pbd_sim_cpp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
