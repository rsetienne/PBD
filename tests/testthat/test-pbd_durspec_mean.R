context("pbd_durspec_mean")

test_that("pbd_durspec_mean and pbd_mean_durspec give same results", {

  # Use parameter names from Etienne & Rosindell 2012
  lambda_3 <- 1.0 # Speciation initiation rate of incipient species
  lambda_2 <-2.0 # Speciation completion rate
  mu_2 <- 0.1 # Extinction rate of incipient species

  durspec_mean_classic <- pbd_durspec_mean(pars = c(lambda_3, lambda_2, mu_2))
  durspec_mean_new <- pbd_mean_durspec(eri = mu_2, scr = lambda_2, siri = lambda_3)
  testthat::expect_equal(durspec_mean_classic, durspec_mean_new)
})
