context("test_PBD")

test_that("Minimal call to PBD", {
  expect_silent(
    pbd_sim(c(0.2, 1, 0.2, 0.1, 0.1), 15)
  )
})

# Moved the pbd_ML test to 'test-pbd_ML.R'
