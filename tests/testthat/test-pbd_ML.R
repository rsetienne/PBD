context("pbd_ML")

test_that("Same likelihood with equivalent conditioning", {
  brts <- 1:10
  initparsopt <- c(4.62,4.34,0.03)

  # pbd_ML will indicate, as warnings, that
  # NAs are produced, as message it will state:
  # Parameter values have been used that cause numerical problems.
  expect_output(
    res1 <- pbd_ML(
        brts = brts, initparsopt = initparsopt,
        exteq = 1,
        cond = 1
    )
  )
  expect_output(
    res2 <- pbd_ML(
        brts = brts, initparsopt = initparsopt,
        exteq = 1,
        cond = 2,
        n_l = 2,
        n_u = Inf
    )
  )
  expect_output(
    res2a <- pbd_ML(
        brts = brts, initparsopt = initparsopt,
        exteq = 1,
        cond = 2,
        n_l = 2,
        n_u = 1000
    )
  )
  expect_equal(res1, res2)
  expect_equal(res2, res2a)
})
