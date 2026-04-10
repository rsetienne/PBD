context("pbd_numspec_vec")

test_that("pbd_numspec_vec2 gives expected results", {

  b_G <- 0.2
  age <- 1
  num <- pbd_numspec_vec2(pars = c(b_G,0,0,0,0),
                          age = age,
                          initvec = c(1, rep(0,92)))
  testthat::expect_equal(num$numspec_tot, 1 + b_G * age, tol = 1E-5)

  b_I <- 0.2
  age <- 1
  num <- pbd_numspec_vec2(pars = c(0,b_I,0,0,0),
                          age = age,
                          initvec = c(rep(0,47),1,rep(0,45)))
  testthat::expect_equal(num$numspec_tot, exp(b_I * age), tol = 1E-5)

  b_I <- 0.2
  mu_I <- 0.1
  age <- 1
  num <- pbd_numspec_vec2(pars = c(0,b_I,0,mu_I,0),
                          age = age,
                          initvec = c(rep(0,47),1,rep(0,45)))
  testthat::expect_equal(num$numspec_tot, exp((b_I - mu_I) * age), tol = 1E-5)

  b_G <- 0.1
  b_I <- b_G
  mu_G <- 0.1
  mu_I <- mu_G
  la <- 0.1
  age <- 1
  num <- pbd_numspec_vec2(pars = c(b_G,b_I,mu_G,mu_I,la),
                          age = age,
                          initvec = c(rep(0,47),1,rep(0,45)))
  testthat::expect_equal(num$numspec_tot, exp((b_G - mu_G) * age), tol = 1E-5)
})
