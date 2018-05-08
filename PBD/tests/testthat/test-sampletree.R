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
    c(1, 2, 3, 1, 1, 1, 0, 0.3, 0.8, -1, 0.9, -1, -1, -1, -1, 1, 2, 1),
    nrow = 3,
    ncol = 6)

  expected_L <- matrix(
    c(3, 2, 1, 1, 1, 1, 0.8, 0.3, 0, -1, 0.9, -1, -1, -1, 1, 1, 2, 1),
    nrow = 3,
    ncol = 6)

  output_L <- PBD:::sampletree(L = input_L, age = 1, samplemethod = "youngest")

  testthat::expect_equal(output_L, expected_L)

})

test_that("sample oldest", {

  input_L <- matrix(
    c(1, 2, 3, 1, 1, 1, 0, 0.3, 0.8, -1, 0.9, -1, -1, -1, -1, 1, 2, 1),
    nrow = 3,
    ncol = 6)

  expected_L <- matrix(
    c(3, 2, 1, 1, 1, 1, 0.8, 0.3, 0, -1, 0.9, -1, 1, -1, -1, 1, 2, 1),
    nrow = 3,
    ncol = 6)

  output_L <- PBD:::sampletree(L = input_L, age = 1, samplemethod = "oldest")

  testthat::expect_equal(output_L, expected_L)

})

test_that("sample shortest", {

  # For students

  # Create an L table that has a different result when sampled
  # with 'shortest' and 'youngest'

  # Create the expected L table

  # Create the L table by sampling 'shortest'

  # Expected and created L table must be identical

})

test_that("sample longest", {

  # For students

  # Create an L table that has a different result when sampled
  # with 'longest' and 'oldest'

  # Create the expected L table

  # Create the L table by sampling 'longest'

  # Expected and created L table must be identical

})

test_that("sample shortest", {

  # For students

  # Create an L table that has a different result when sampled
  # with 'shortest' and 'youngest'

  # Create the expected L table

  # Create the L table by sampling 'shortest'

  # Expected and created L table must be identical

})

test_that("abuse", {

  # For students

  # We cannot change the interface of 'sampletree',
  # i.e. we cannot add if-then-stop statements there

  # We can explicitly test for the current errors
  # 'expect_error' with a negative crown age, what error does it give now?
  # Test for those errors explicity.
})
