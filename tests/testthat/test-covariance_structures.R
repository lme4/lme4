library(lme4) 
library(testthat)

# UNIT TESTS FOR DiagonalCov

context("Testing DiagonalCov Class and Methods")

test_that("DiagonalCov object creation and initialization works", {
  dim <- 3
  diag_cov_obj <- new("DiagonalCov", dimension = dim)

  expect_s4_class(diag_cov_obj, "DiagonalCov")
  expect_equal(diag_cov_obj@dimension, dim)

  start_vals <- get_start_values(diag_cov_obj)
  expect_equal(length(start_vals), dim)
  expect_equal(start_vals, rep(log(1), dim))
})

test_that("DiagonalCov parameter setting and getting works", {
  dim <- 4
  diag_cov_obj <- new("DiagonalCov", dimension = dim)
  
  # Test set_parameters
  unconstrained_params <- log(c(1.5, 2.0, 0.5, 3.0))
  diag_cov_obj <- set_parameters(diag_cov_obj, unconstrained_params)
  
  expect_equal(diag_cov_obj@internal_diag_values, c(1.5, 2.0, 0.5, 3.0))
  expect_equal(diag_cov_obj@parameters, unconstrained_params) 
  expect_error(set_parameters(diag_cov_obj, c(1, 2)), 
               "Incorrect number of parameters for DiagonalCov. Expected 4")
  
  # Test get_parameters
  expect_equal(get_parameters(diag_cov_obj), unconstrained_params)
})


test_that("DiagonalCov matrix computations are correct", {
    dim <- 3
    diag_cov_obj <- new("DiagonalCov", dimension = dim)
    variances <- c(2, 4, 5)
    diag_cov_obj <- set_parameters(diag_cov_obj, log(variances))

    # compute_covariance_matrix
    expected_cov <- diag(variances)
    computed_cov <- compute_covariance_matrix(diag_cov_obj)
    expect_equal(computed_cov, expected_cov)

    # compute_log_det_covariance_matrix
    expected_log_det <- sum(log(variances))
    computed_log_det <- compute_log_det_covariance_matrix(diag_cov_obj)
    expect_equal(computed_log_det, expected_log_det)

    # compute_inverse_covariance_matrix
    expected_inv_cov <- diag(1 / variances)
    computed_inv_cov <- compute_inverse_covariance_matrix(diag_cov_obj)
    expect_equal(computed_inv_cov, expected_inv_cov)

    # get_cholesky_factor
    expected_chol <- diag(sqrt(variances))
    computed_chol <- get_cholesky_factor(diag_cov_obj)
    expect_equal(computed_chol, expected_chol)
})

test_that("Diagonal validation works", {
    dim <- 3
    diag_cov_obj <- new("DiagonalCov", dimension = dim)

    # Validity check
    valid_params <- log(c(1, 2, 3))
    diag_cov_obj <- set_parameters(diag_cov_obj, valid_params)
    expect_true(validate_parameters(diag_cov_obj))

    # Invalidity Check 
    invalid_params <- log(c(1, 0, 3))
    diag_cov_obj_invalid <- set_parameters(diag_cov_obj, invalid_params)
    expect_error(validate_parameters(diag_cov_obj_invalid),
                 "DiagonalCov: All variances must be positive.")
})

test_that("DiagonalCov show method executes without error", {
    dim <- 2
    diag_cov_obj <- new("DiagonalCov", dimension = dim)
    diag_cov_obj <- set_parameters(diag_cov_obj, log(c(10, 20)))

    expect_output(show(diag_cov_obj), "DiagonalCov object \\(dimension: 2\\)")
    expect_output(show(diag_cov_obj), "Variances: 10 20")
})


# EDGE CASE TESTS FOR DiagonalCov

context("Testing DiagonalCov Edge Cases")

test_that("DiagonalCov works correctly for dimension = 1", {
    dim <- 1
    diag_cov_obj <- new("DiagonalCov", dimension = dim)
    diag_cov_obj <- set_parameters(diag_cov_obj, log(4.0))

    expect_equal(n_parameters(diag_cov_obj), 1)
    expect_equal(compute_covariance_matrix(diag_cov_obj), matrix(4.0))
    expect_equal(compute_log_det_covariance_matrix(diag_cov_obj), log(4.0))
    expect_equal(get_cholesky_factor(diag_cov_obj), matrix(2.0))
})

test_that("DiagonalCov handles dimension = 0 in a sensible way", {
    dim <- 0
    diag_cov_obj <- new("DiagonalCov", dimension = dim)

    expect_equal(n_parameters(diag_cov_obj), 0)
    expect_equal(get_start_values(diag_cov_obj), numeric(0))
    
    # The covariance matrix should be a 0x0 matrix
    cov_mat <- compute_covariance_matrix(diag_cov_obj)
    expect_equal(nrow(cov_mat), 0)
    expect_equal(ncol(cov_mat), 0)

    # The log-determinant of a 0-dim matrix should be 0
    expect_equal(compute_log_det_covariance_matrix(diag_cov_obj), 0)
})

test_that("DiagonalCov methods handle non-finite parameters", {
    dim <- 3
    diag_cov_obj <- new("DiagonalCov", dimension = dim)

    # Test with Inf (valid)
    params_inf <- c(log(2), Inf, log(4))
    diag_cov_obj_inf <- set_parameters(diag_cov_obj, params_inf)
    expect_true(validate_parameters(diag_cov_obj_inf)) # Inf > 0 is TRUE
    expect_equal(compute_log_det_covariance_matrix(diag_cov_obj_inf), Inf)
  
    # Test with NA (invalid)
    params_na <- c(log(2), NA, log(4))
    diag_cov_obj_na <- set_parameters(diag_cov_obj, params_na)
    expect_error(validate_parameters(diag_cov_obj_na))
})
