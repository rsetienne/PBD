## ------------------------------------------------------------------------
library(PBD)

## ------------------------------------------------------------------------
library(ape)

## ------------------------------------------------------------------------
seed <- 43
set.seed(seed)
b_1 <- 0.3 # speciation-initiation rate of good species
la_1 <- 0.2 # speciation-completion rate
b_2 <- b_1 # the speciation-initiation rate of incipient species
mu_1 <- 0.1 #  extinction rate of good species
mu_2 <- mu_1 # extinction rate of incipient species 
pars <- c(b_1, la_1, b_2, mu_1, mu_2)
age <- 15 # the age for the simulation 
pbd_sim_result <- pbd_sim(pars = pars, age = age, plotit = TRUE)
phylogeny <- pbd_sim_result$recontree
plot(phylogeny)
plot(pbd_sim_result$igtree.extant)
plot(pbd_sim_result$tree)
names(pbd_sim_result)

## ------------------------------------------------------------------------
n_taxa <- length(phylogeny$tip.label)
print(paste("The phylogeny has", n_taxa, "taxa"))

## ------------------------------------------------------------------------
brts <- branching.times(phylogeny)  # branching times
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

reltolx <- 10^-6 # relative tolerance of parameter values in optimization
reltolf <- 10^-6 # relative tolerance of function value in optimization
abstolx <- 10^-6 # absolute tolerance of parameter values in optimization 
tol <- c(reltolx, reltolf, abstolx)

## ------------------------------------------------------------------------
r <- pbd_ML(
  brts = brts,
  initparsopt = initparsopt, 
  exteq = exteq,
  soc = soc, 
  cond = cond,
  btorph = btorph,
  tol = tol,
  verbose = FALSE
)

## ------------------------------------------------------------------------
print(r)

## ------------------------------------------------------------------------
df <- as.data.frame(x = list(b = b_1, mu_1 = mu_1, lambda_1 = la_1, mu_2 = mu_2,  loglik = NA, df = NA, conv = NA))
df <- rbind(df, r)
row.names(df) <- c("true", "ML")
knitr::kable(df)

