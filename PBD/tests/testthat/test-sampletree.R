context("sampletree")

test_that("sample randomly", {

  # For students

  # Hints/requirements:
  # - Use PBD:::sampletree
  # - lintr-bot prefers l_table over L

  # Create an L table that can be sampled randomly in two ways

  # Create the two expected L tables

  # 10x: create a randomly sampled L table
  #   - each of these must be one of the two expected L tables

})

test_that("sample youngest", {

  input_L <- matrix(
    c(1, 2, 3, 4, 0, 1, 2, 1, -1e-10, 0.2, 0.4, 0.6, 0, 1.0, -1, 0.8, -1, -1, -1, -1, 1, 2, 1, 3),
    nrow = 4,
    ncol = 6)

  expected_L <- matrix(
    c(4, 3, 2, 1, 1, 2, 1, 0, 0.6, 0.4, 0.2, -1e-10, 0.8, -1, 1, 0, -1, -1, -1, 1.2, 3, 1, 2, 1),
    nrow = 4,
    ncol = 6)

  output_L <- PBD:::sampletree(L = input_L, age = 1.2, samplemethod = "youngest")

  testthat::expect_equal(output_L, expected_L)

})

test_that("sample oldest", {

  input_L <- matrix(
    c(1, 2, 3, 4, 0, 1, 2, 1, -1e-10, 0.2, 0.4, 0.6, 0, 1.0, -1, 0.8, -1, -1, -1, -1, 1, 2, 1, 3),
    nrow = 4,
    ncol = 6)

  expected_L <- matrix(
    c(4, 3, 2, 1, 1, 2, 1, 0, 0.6, 0.4, 0.2, -1e-10, 0.8, -1, 1, 0, -1, 1.2, -1, -1, 3, 1, 2, 1),
    nrow = 4,
    ncol = 6)

  output_L <- PBD:::sampletree(L = input_L, age = 1.2, samplemethod = "oldest")

  testthat::expect_equal(output_L, expected_L)

})

test_that("sample shortest", {

  #           +------------------------- 1-3
  #           |
  #     +-----+
  #     |     |
  # ====+     +-----------------=======  2-2
  #     |
  #     |            +----=============  3-4
  #     |            |
  #     +============+
  #                  |
  #                  +=================  1-1
  #
  # +----+-----+-----+-----+----+-----+- t
  # 0   0.2   0.4   0.6   0.8  1.0   1.2
  #
  input_L <- matrix(
    c(
      1, 2, 3, 4, # Subspecies IDs
      0, 1, 2, 1, # Parent IDs
      -1e-10, 0.2, 0.4, 0.6, # Speciation initiation times
      0, 1.0, -1, 0.8, # Speciation completion times, -1 if still incipient
      -1, -1, -1, -1, # Extinction time, -1 if still present
      1, 2, 1, 3 # Species ID
    ),
    nrow = 4,
    ncol = 6)

  expected_L <- matrix(
    c(4, 3, 2, 1, 1, 2, 1, 0, 0.6, 0.4, 0.2, -1e-10, 0.8, -1, 1, 0, -1, 1.2, -1, -1, 3, 1, 2, 1),
    nrow = 4,
    ncol = 6)

  output_L <- PBD:::sampletree(L = input_L, age = 1.2, samplemethod = "shortest")

  testthat::expect_equal(output_L, expected_L)

})

test_that("sample longest", {

  input_L <- matrix(
    c(1, 2, 3, 4, 0, 1, 2, 1, -1e-10, 0.2, 0.4, 0.6, 0, 1.0, -1, 0.8, -1, -1, -1, -1, 1, 2, 1, 3),
    nrow = 4,
    ncol = 6)

  expected_L <- matrix(
    c(4, 3, 2, 1, 1, 2, 1, 0, 0.6, 0.4, 0.2, -1e-10, 0.8, -1, 1, 0, -1, -1, -1, 1.2, 3, 1, 2, 1),
    nrow = 4,
    ncol = 6)

  output_L <- PBD:::sampletree(L = input_L, age = 1.2, samplemethod = "longest")

  testthat::expect_equal(output_L, expected_L)

})

test_that("abuse", {

  # For students

  # We cannot change the interface of 'sampletree',
  # i.e. we cannot add if-then-stop statements there

  # We can explicitly test for the current errors
  # 'expect_error' with a negative crown age, what error does it give now?
  # Test for those errors explicity.
})
