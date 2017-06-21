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

test_that("siri (lambda_3) is zero, equation 21a", {

  # Equation 21a in
  #' Etienne, Rampal S., and James Rosindell. "Prolonging the past
  #'   counteracts the pull of the present: protracted speciation can explain
  #'   observed slowdowns in diversification." Systematic
  #'   Biology 61.2 (2012): 204-213.
  eri <- 0.1 # extinction rate of incipient species
  scr <- 0.2 # speciation completion rate
  siri <- 0.0 # speciation initiation rate of incipient species
  mean_durspec <- PBD::pbd_mean_durspec(eri = eri, scr = scr, siri = siri)

  # Calculate by hand
  mu_2 <- eri
  lambda_2 <- scr
  lambda_3 <- siri
  expected_mean_durspec <- 1.0 / (lambda_2 + mu_2)

  expected_mean_durspec_too <- 3.333333
  testthat::expect_equal(mean_durspec, expected_mean_durspec,
    tolerance = 0.000001)
  testthat::expect_equal(mean_durspec, expected_mean_durspec_too,
    tolerance = 0.000001)
})

test_that("eri (mu_2) is zero, equation 21b", {

  # Equation 21b in
  #' Etienne, Rampal S., and James Rosindell. "Prolonging the past
  #'   counteracts the pull of the present: protracted speciation can explain
  #'   observed slowdowns in diversification." Systematic
  #'   Biology 61.2 (2012): 204-213.
  eri <- 0.0 # extinction rate of incipient species
  scr <- 0.2 # speciation completion rate
  siri <- 0.3 # speciation initiation rate of incipient species
  mean_durspec <- PBD::pbd_mean_durspec(eri = eri, scr = scr, siri = siri)

  # Calculate by hand
  mu_2 <- eri
  lambda_2 <- scr
  lambda_3 <- siri
  expected_mean_durspec <- (1.0 / lambda_3) * log(1.0 + (lambda_3 / lambda_2))

  expected_mean_durspec_too <- 3.054302
  testthat::expect_equal(mean_durspec, expected_mean_durspec,
    tolerance = 0.000001)
  testthat::expect_equal(mean_durspec, expected_mean_durspec_too,
    tolerance = 0.000001)
})

test_that("complete calculation", {

  eri <- 0.2  # extinction rate of incipient species, mu_2
  scr <- 0.3  # speciation completion rate, lambda_2
  siri <- 0.5 # speciation initiation rate of incipient species, lambda_3

  mean_durspec <- PBD::pbd_mean_durspec(eri = eri, scr = scr, siri = siri)

  # Calculate by hand
  mu_2 <- eri
  lambda_2 <- scr
  lambda_3 <- siri
  d_term_1 <- (lambda_2 + lambda_3) ^ 2
  d_term_2 <- 2.0 * (lambda_2 - lambda_3) * mu_2
  d_term_3 <- mu_2 ^ 2
  d <- sqrt(d_term_1 + d_term_2 + d_term_3)
  tau_term_1 <- 2.0 / (d - lambda_2 + lambda_3 - mu_2)
  tau_term_2_fraction <- (lambda_2 - lambda_3 + mu_2) / d
  tau_term_2 <- 2.0 / (1.0 + tau_term_2_fraction)
  expected_mean_durspec <- tau_term_1 * log(tau_term_2)

  expected_mean_durspec_too <- 1.789698

  testthat::expect_equal(mean_durspec, expected_mean_durspec,
    tolerance = 0.000001)
  testthat::expect_equal(mean_durspec, expected_mean_durspec_too,
    tolerance = 0.000001)
})


test_that("abuse", {

  eri <- 0.2  # extinction rate of incipient species, mu_2
  scr <- 0.3  # speciation completion rate, lambda_2
  siri <- 0.5 # speciation initiation rate of incipient species, lambda_3

  testthat::expect_error(
    PBD::pbd_mean_durspec(eri = -1.2, scr = scr, siri = siri),
    "extinction rate of incipient species must be zero or positive"
  )
  testthat::expect_error(
    PBD::pbd_mean_durspec(eri = NaN, scr = scr, siri = siri),
    "extinction rate of incipient species must be zero or positive"
  )

  testthat::expect_error(
    PBD::pbd_mean_durspec(eri = eri, scr = -1.2, siri = siri),
    "speciation completion rate must be zero or positive"
  )
  testthat::expect_error(
    PBD::pbd_mean_durspec(eri = eri, scr = NaN, siri = siri),
    "speciation completion rate must be zero or positive"
  )

  testthat::expect_error(
    PBD::pbd_mean_durspec(eri = eri, scr = scr, siri = -1.2),
    "speciation initiation rate of incipient species must be zero or positive"
  )
  testthat::expect_error(
    PBD::pbd_mean_durspec(eri = eri, scr = scr, siri = NaN),
    "speciation initiation rate of incipient species must be zero or positive"
  )
})

test_that("basic vectorized usage", {

  eri_1 <- 0.1 # incipient species extinction rate
  scr_1 <- 0.2 # speciation completion rate
  siri_1 <- 0.3 # speciation initiation rate of incipient species
  eri_2 <- 0.4 # incipient species extinction rate
  scr_2 <- 0.5 # speciation completion rate
  siri_2 <- 0.6 # speciation initiation rate of incipient species

  mean_durspec_1 <- PBD::pbd_mean_durspec(eri = eri_1, scr = scr_1, siri = siri_1)
  mean_durspec_2 <- PBD::pbd_mean_durspec(eri = eri_2, scr = scr_2, siri = siri_2)
  mean_durspecs_expected <- c(mean_durspec_1, mean_durspec_2)

  eris <- c(eri_1, eri_2)
  scrs <- c(scr_1, scr_2)
  siris <- c(siri_1, siri_2)

  mean_durspecs <- PBD::pbd_mean_durspecs(eri = eris, scr = scrs, siri = siris)

  testthat::expect_equal(mean_durspecs, mean_durspecs_expected,
    tolerance = 0.000001)
})

test_that("vectorized usage with recycling", {

  # Recycle eri
  testthat::expect_silent(
    PBD::pbd_mean_durspecs(eri = c(0.1), scr = c(0.2, 0.3), siri = c(0.4, 0.5))
  )
  # Recycle scr
  testthat::expect_silent(
    PBD::pbd_mean_durspecs(eri = c(0.1, 0.2), scr = c(0.3), siri = c(0.4, 0.5))
  )
  # Recycle siri
  testthat::expect_silent(
    PBD::pbd_mean_durspecs(eri = c(0.1, 0.2), scr = c(0.3, 0.4), siri = c(0.5))
  )
})

test_that("vectorized usage with invalid arguments", {

  testthat::expect_error(
    PBD::pbd_mean_durspecs(
      eri = c(-0.1, 0.2), # Invalid first eri
      scr = c(0.3, 0.4),
      siri = c(0.5, 0.6)
    )
  )

  testthat::expect_error(
    PBD::pbd_mean_durspecs(
      eri = c(0.1, 0.2),
      scr = c(0.3, -0.4), # Invalid second scr
      siri = c(0.5, 0.6)
    )
  )

  testthat::expect_error(
    PBD::pbd_mean_durspecs(
      eri = c(-0.1, 0.2),
      scr = c(0.3, 0.4),
      siri = c(NA, 0.6) # Invalid first siri
    )
  )

})

