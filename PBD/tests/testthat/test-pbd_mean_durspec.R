context("pbd_mean_durspec")

test_that("example", {

  eri <- 0.1 # incipient species extinction rate
  scr <- 0.2 # speciation completion rate
  siri <- 0.3 # speciation initiation rate of incipient species
  mean_durspec <- PBD::pbd_mean_durspec(eri = eri, scr = scr, siri = siri)
  expected_mean_durspec <- 2.829762
  testthat::expect_equal(mean_durspec, expected_mean_durspec,
    tolerance = 0.000001)

})
