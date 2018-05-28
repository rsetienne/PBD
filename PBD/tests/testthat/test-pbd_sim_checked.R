context("pbd_sim_checked")

test_that("use", {

  expect_silent(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = 1.0
    )
  )
})

test_that("PBD::pbd_sim and pbd_sim_checked must give same results", {

  set.seed(42)
  checked <- pbd_sim_checked(
    erg = 0.0,
    eri = 0.0,
    scr = 1.0,
    sirg = 1.0,
    siri = 1.0,
    crown_age = 1.0
  )

  set.seed(42)
  normal <- PBD::pbd_sim(
    pars = c(1.0, 1.0, 1.0, 0.0, 0.0),
    age = 1.0,
    soc = 2
  )

  expect_equal(checked, normal)

})

test_that("abuse", {

  expect_error(
    pbd_sim_checked(
      erg = -123, # Error
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = 1.0
    ),
    "'erg' must be positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = -123, # Error
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = 1.0
    ),
    "'eri' must be positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = -123, # Error
      sirg = 1.0,
      siri = 1.0,
      crown_age = 1.0
    ),
    "'scr' must be positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = -123, # Error
      siri = 1.0,
      crown_age = 1.0
    ),
    "'sirg' must be positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = -123, #Error
      crown_age = 1.0
    ),
    "'siri' must be positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = NULL,
      stem_age = NULL
    ),
    paste0(
      "At least one of 'crown_age' or 'stem_age' ",
      "must be non-zero and positive"
    )
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = -123,
      stem_age = NULL
    ),
    "'crown_age' must be non-zero and positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = NULL,
      stem_age = -123
    ),
    "'stem_age' must be non-zero and positive"
  )

  expect_error(
    pbd_sim_checked(
      erg = 0.0,
      eri = 0.0,
      scr = 1.0,
      sirg = 1.0,
      siri = 1.0,
      crown_age = 123,
      stem_age = 234
    ),
    "Must set either 'crown_age' or 'stem_age'"
  )

})
