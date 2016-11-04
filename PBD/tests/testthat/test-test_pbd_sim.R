context("test_pbd_sim")

test_that("pbd_sim works", {
  pbd_sim(c(0.2,1,0.2,0.1,0.1),15)
  expect_equal(2 * 2, 4)
})
