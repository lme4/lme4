library(lme4) 
library(testthat)

# Tests for DiagonalCov

context("Testing DiagonalCov Class and Methods")

test_that("DiagonalCov object creation and initialization works", {
	dim <- 3L
	params <- rep(log(1), dim)
	diag_cov_obj <- new("DiagonalCov", dimension = dim, parameters = params)

	expect_s4_class(diag_cov_obj, "DiagonalCov")
	expect_true(validObject(diag_cov_obj))
	expect_equal(diag_cov_obj@dimension, dim)
	expect_equal(get_parameters(diag_cov_obj), params)

  	# Test for invalid creation
	expect_error(new("DiagonalCov", dimension = dim, parameters = rep(0, dim - 1)))
})

test_that("DiagonalCov parameter setting and metadata methods work", {
	dim <- 3L
  	params <- rep(log(1), dim)
  	obj <- new("DiagonalCov", dimension = dim, parameters = params)

  	new_params <- c(1, 2, 3)
  	obj <- set_parameters(obj, new_params)
  	expect_equal(get_parameters(obj), new_params)
  	expect_error(set_parameters(obj, c(1, 2))) # Incorrect length

  	expect_equal(get_start_values(obj), rep(log(1), dim))
  	expect_true(is_diagonal(obj))
  	expect_equal(get_lower_bounds(obj), rep(-Inf, dim))
})

test_that("DiagonalCov computations are correct", {
	dim <- 3L
	log_variances <- log(c(2, 4, 5))
	diag_cov_obj <- new("DiagonalCov", dimension = dim, parameters = log_variances)

  	interp_params <- get_interpretable_parameters(diag_cov_obj)
  	expect_equal(interp_params$variances, c(2, 4, 5))

  	expect_equal(as.matrix(compute_covariance_matrix(diag_cov_obj)), diag(c(2, 4, 5)))
  	expect_equal(as.matrix(get_cholesky_factor(diag_cov_obj)), diag(sqrt(c(2, 4, 5))))
  	expect_equal(compute_log_det_covariance_matrix(diag_cov_obj), sum(log_variances))
})

context("Testing DiagonalCov Edge Cases")

test_that("DiagonalCov works correctly for dimension = 1", {
  	d <- 1L
  	obj <- new("DiagonalCov", dimension = d, parameters = log(4))
  	expect_equal(as.matrix(compute_covariance_matrix(obj)), matrix(4))
  	expect_equal(compute_log_det_covariance_matrix(obj), log(4))
})

test_that("DiagonalCov handles dimension = 0 gracefully", {
  	d <- 0L
  	obj <- new("DiagonalCov", dimension = d, parameters = numeric(0))
  	expect_true(validObject(obj))
  	expect_equal(nrow(compute_covariance_matrix(obj)), 0)
  	expect_equal(compute_log_det_covariance_matrix(obj), 0)
})

# Tests for UnstructuredCov

context("Testing UnstructuredCov Class and Methods")

test_that("UnstructuredCov object creation and validity works", {
  	d <- 3L
  	n_params <- d * (d + 1) / 2
  	params <- rep(0, n_params)
  	obj <- new("UnstructuredCov", dimension = d, parameters = params)

  	expect_s4_class(obj, "UnstructuredCov")
  	expect_true(validObject(obj))

  	expect_error(new("UnstructuredCov", dimension = d, parameters = rep(0, n_params - 1)))
})

test_that("UnstructuredCov parameter getting/setting works", {
	d <- 2L
  	n_params <- d * (d + 1) / 2
  	obj <- new("UnstructuredCov", dimension = d, parameters = rep(0, n_params))

  	true_params <- c(log(2), 0.5, log(3)) 
  	obj <- set_parameters(obj, true_params)
  	expect_equal(get_parameters(obj), true_params)

  	L_matrix <- matrix(c(2, 0.5, 0, 3), nrow = 2) 
  	obj_from_mat <- set_parameters(obj, L_matrix)
  	expect_equal(get_parameters(obj_from_mat), true_params)

  	expect_error(set_parameters(obj, c(1, 2))) # wrong length
  	expect_error(set_parameters(obj, matrix(1:9, 3, 3))) # wrong dimension
})

test_that("UnstructuredCov matrix computations are correct", {
	d <- 2L
	# Corresponds to L = matrix(c(2, 0.5, 0, 3), 2)
  
	params <- c(log(2), 0.5, log(3))
	obj <- new("UnstructuredCov", dimension = d, parameters = params)

	L_expected <- matrix(c(2, 0.5, 0, 3), nrow = 2)
	Sigma <- L_expected %*% t(L_expected)

   	expect_equal(as.matrix(get_cholesky_factor(obj)), L_expected)

    	expect_equal(as.matrix(compute_covariance_matrix(obj)), Sigma)

      	expect_equal(compute_log_det_covariance_matrix(obj), 2 * (log(2) + log(3)))

    	expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), solve(Sigma))
})

test_that("UnstructuredCov metadata and show method work", {
  	d <- 3L
  	n_params <- d * (d + 1) / 2
  	obj <- new("UnstructuredCov", dimension = d, parameters = rep(0, n_params))

  	expect_equal(n_parameters(obj), n_params)
  	expect_false(is_diagonal(obj))
  	expect_equal(get_lower_bounds(obj), rep(-Inf, n_params))

   	params <- get_interpretable_parameters(obj)
  	expect_equal(params$st_devs, c(1, 1, 1))
  	expect_equal(as.matrix(params$covariance_matrix), diag(3))

    	expect_output(show(obj), "UnstructuredCov object \\(dimension: 3\\)")
  	expect_output(show(obj), "Standard Deviations: 1 1 1")
})


context("Testing UnstructuredCov Edge Cases")

test_that("UnstructuredCov works correctly for dimension = 1", {
  	d <- 1L
  	log_L11 <- log(2)
  	obj <- new("UnstructuredCov", dimension = d, parameters = log_L11)
  	expect_true(validObject(obj))

  	expected_variance <- exp(log_L11)^2
  	expect_equal(as.matrix(compute_covariance_matrix(obj)), matrix(expected_variance))

  	expected_log_det <- 2 * log_L11
  	expect_equal(compute_log_det_covariance_matrix(obj), expected_log_det)
})

test_that("UnstructuredCov handles dimension = 0", {
  	d <- 0L
  	obj <- new("UnstructuredCov", dimension = d, parameters = numeric(0))
  	expect_true(validObject(obj))
  	expect_equal(nrow(get_cholesky_factor(obj)), 0)
  	expect_equal(nrow(compute_covariance_matrix(obj)), 0)
  	expect_equal(compute_log_det_covariance_matrix(obj), 0)
})

test_that("UnstructuredCov handles non-finite parameters", {
  	d <- 2L
  	# NA in off-diagonal
  	params_na <- c(log(2), NA, log(3))
  	obj_na <- new("UnstructuredCov", dimension = d, parameters = params_na)
  	expect_true(is.na(compute_covariance_matrix(obj_na)[1, 2]))

  	# Inf in off-diagonal
  	params_inf <- c(log(2), Inf, log(3))
  	obj_inf <- new("UnstructuredCov", dimension = d, parameters = params_inf)
  	expect_true(is.infinite(compute_covariance_matrix(obj_inf)[1, 2]))
})

test_that("UnstructuredCov handles singular matrices", {
  	d <- 2L
  	params_singular <- c(-Inf, 0.5, log(3))
  	obj_singular <- new("UnstructuredCov", dimension = d, parameters = params_singular)

  	expect_equal(compute_log_det_covariance_matrix(obj_singular), -Inf)

  	# Inverse should not be computable
  	expect_error(compute_inverse_covariance_matrix(obj_singular), "singular")
})

