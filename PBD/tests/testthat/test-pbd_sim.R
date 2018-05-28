context("pbd_sim")

test_that("Default run gives default output", {

  out <- pbd_sim(
    pars = c(2, 0.1, 2, 0, 0),
    age = 1
  )

  testthat::expect_equal(length(out), 13)
  expected_names <- c(
    "tree", "stree_random", "stree_oldest", "stree_youngest", "L",
    "sL_random", "sL_oldest", "sL_youngest", "igtree.extinct",
    "igtree.extant", "recontree", "reconL", "L0"
  )
  testthat::expect_equal(names(out), expected_names)
})

test_that("Run with shortest/longest gives shortest/longest output", {

  out <- pbd_sim(
    pars = c(2, 0.1, 2, 0, 0),
    age = 1,
    add_shortest_and_longest = TRUE)

  testthat::expect_equal(length(out), 17)
  expected_names <- c(
    "tree", "stree_random", "stree_oldest", "stree_youngest", "L", "sL_random",
    "sL_oldest", "sL_youngest", "igtree.extinct", "igtree.extant", "recontree",
    "reconL", "L0", "stree_shortest", "stree_longest", "sL_shortest",
    "sL_longest"
  )
  testthat::expect_equal(names(out), expected_names)
})
