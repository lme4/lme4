library(testthat)
library(Matrix)

# Helper to allow tests to extract coordinates from a simple list
setMethod("get_coordinates", "list", function(data_context) data_context$coords)

# Tests for FlexibleCovariance with Diagonal Structure

context("Testing FlexibleCovariance with Diagonal Structure")

test_that("Diagonal model object creation and initialization works", {
	d <- 3L
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = log(c(1, 2, 3)))
	cor_model <- new("IdentityCorrelationModel", parameters = numeric(0))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	expect_s4_class(obj, "FlexibleCovariance")
	expect_true(validObject(obj))
	expect_equal(obj@dimension, d)
	expect_equal(get_parameters(obj), log(c(1, 2, 3)))

	# Test for invalid creation (error is in the component model)
	expect_error(
		new("HeterogeneousVarianceModel", dimension = d, parameters = log(c(1, 2))),
		"Number of parameters must equal dimension."
	)
})

test_that("Diagonal model parameter setting and metadata methods work", {
	d <- 3L
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = log(c(1, 2, 3)))
	cor_model <- new("IdentityCorrelationModel", parameters = numeric(0))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	new_params <- c(1, 2, 3)
	obj <- set_parameters(obj, new_params)
	expect_equal(get_parameters(obj), new_params)
	expect_error(set_parameters(obj, c(1, 2))) # Incorrect length

	expect_equal(get_start_values(obj), rep(0, d))
	expect_true(is_diagonal(obj))
	expect_equal(get_lower_bounds(obj), rep(-Inf, d))
})

test_that("Diagonal model computations are correct", {
	d <- 3L
	st_devs <- c(2, 3, 4)
	params <- log(st_devs^2)
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = params)
	cor_model <- new("IdentityCorrelationModel", parameters = numeric(0))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	interp_params <- get_interpretable_parameters(obj)
	expect_equal(interp_params$st_devs, st_devs)

	Sigma_expected <- Diagonal(d, x = st_devs^2)
	expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma_expected))
	expect_equal(as.matrix(get_cholesky_factor(obj)), as.matrix(diag(st_devs)))
	expect_equal(compute_log_det_covariance_matrix(obj), sum(params))
	expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), as.matrix(diag(1/st_devs^2)))
})

test_that("show method for Diagonal model works", {
	d <- 3L
	st_devs <- c(2, 3.5, 4)
	params <- log(st_devs^2)
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = params)
	cor_model <- new("IdentityCorrelationModel", parameters = numeric(0))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	expect_output(show(obj), "Heterogeneous-Identity Covariance object \\(dimension: 3\\)")
	expect_output(show(obj), "  Standard Deviations: 2\\.0000 3\\.5000 4\\.0000")
})

context("Testing FlexibleCovariance (Diagonal) Edge Cases")

test_that("Diagonal model works correctly for dimension = 1", {
	d <- 1L
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = log(4))
	cor_model <- new("IdentityCorrelationModel")
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	expect_equal(as.matrix(compute_covariance_matrix(obj)), matrix(4))
	expect_equal(compute_log_det_covariance_matrix(obj), log(4))
})

test_that("Diagonal model handles dimension = 0", {
	d <- 0L
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = numeric(0))
	cor_model <- new("IdentityCorrelationModel")
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	expect_true(validObject(obj))
	expect_equal(nrow(compute_covariance_matrix(obj)), 0)
	expect_equal(compute_log_det_covariance_matrix(obj), 0)
})

# Tests for UnstructuredCovariance

context("Testing UnstructuredCovariance Class and Methods")

test_that("UnstructuredCovariance object creation and validity works", {
	d <- 3L
	n_params <- d * (d + 1) / 2
	params <- rep(0, n_params)
	obj <- new("UnstructuredCovariance", dimension = d, parameters = params)

	expect_s4_class(obj, "UnstructuredCovariance")
	expect_true(validObject(obj))

	expect_error(new("UnstructuredCovariance", dimension = d, parameters = rep(0, n_params - 1)))
})

test_that("UnstructuredCovariance parameter getting/setting works", {
	d <- 2L
	n_params <- d * (d + 1) / 2
	obj <- new("UnstructuredCovariance", dimension = d, parameters = rep(0, n_params))

	true_params <- c(log(2), 0.5, log(3))
	obj <- set_parameters(obj, true_params)
	expect_equal(get_parameters(obj), true_params)

	expect_error(set_parameters(obj, c(1, 2))) # wrong length
})

test_that("UnstructuredCovariance matrix computations are correct", {
	d <- 2L
	# Corresponds to L = matrix(c(2, 0.5, 0, 3), 2)
	params <- c(log(2), 0.5, log(3))
	obj <- new("UnstructuredCovariance", dimension = d, parameters = params)

	L_expected <- matrix(c(2, 0.5, 0, 3), nrow = 2)
	Sigma <- L_expected %*% t(L_expected)

	expect_equal(as.matrix(get_cholesky_factor(obj)), L_expected)
	expect_equal(as.matrix(compute_covariance_matrix(obj)), Sigma)
	expect_equal(compute_log_det_covariance_matrix(obj), 2 * (log(2) + log(3)))
	expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), solve(Sigma))
})

test_that("UnstructuredCovariance metadata methods work", {
	d <- 3L
	n_params <- d * (d + 1) / 2
	obj <- new("UnstructuredCovariance", dimension = d, parameters = rep(0, n_params))

	expect_equal(n_parameters(obj), n_params)
	expect_false(is_diagonal(obj))
	expect_equal(get_lower_bounds(obj), rep(-Inf, n_params))

	params <- get_interpretable_parameters(obj)
	expect_equal(params$st_devs, c(1, 1, 1))
})

test_that("show method for UnstructuredCovariance works", {
	d <- 2L
	# Corresponds to L = matrix(c(2, 0.5, 0, 3), 2)
	params <- c(log(2), 0.5, log(3))
	obj <- new("UnstructuredCovariance", dimension = d, parameters = params)

	expect_output(show(obj), "UnstructuredCovariance object \\(dimension: 2\\)")
	expect_output(show(obj), "  Standard Deviations: 2\\.0000 3\\.0414")
})

context("Testing UnstructuredCovariance Edge Cases")

test_that("UnstructuredCovariance works correctly for dimension = 1", {
	d <- 1L
	log_L11 <- log(2)
	obj <- new("UnstructuredCovariance", dimension = d, parameters = log_L11)
	expect_true(validObject(obj))

	expected_variance <- exp(log_L11)^2
	expect_equal(as.matrix(compute_covariance_matrix(obj)), matrix(expected_variance))

	expected_log_det <- 2 * log_L11
	expect_equal(compute_log_det_covariance_matrix(obj), expected_log_det)
})

test_that("UnstructuredCovariance handles dimension = 0", {
	d <- 0L
	obj <- new("UnstructuredCovariance", dimension = d, parameters = numeric(0))
	expect_true(validObject(obj))
	expect_equal(nrow(get_cholesky_factor(obj)), 0)
	expect_equal(nrow(compute_covariance_matrix(obj)), 0)
	expect_equal(compute_log_det_covariance_matrix(obj), 0)
})

test_that("UnstructuredCovariance handles singular matrices", {
	d <- 2L
	params_singular <- c(-Inf, 0.5, log(3))
	obj_singular <- new("UnstructuredCovariance", dimension = d, parameters = params_singular)

	expect_equal(compute_log_det_covariance_matrix(obj_singular), -Inf)
	expect_error(compute_inverse_covariance_matrix(obj_singular), "singular")
})

# Tests for FlexibleCovariance with Homogeneous CS Structure

context("Testing FlexibleCovariance with Homogeneous CS Structure")

test_that("Homogeneous-CS model object creation and validity works", {
	d <- 3L
	st_dev <- 2
	rho <- 0.5
	var_model <- new("HomogeneousVarianceModel", parameters = log(st_dev^2))
	cor_model <- new("CSCorrelationModel", parameters = atanh(rho))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	expect_s4_class(obj, "FlexibleCovariance")
	expect_true(validObject(obj))
	expect_equal(obj@dimension, d)

	# Test invalid rho
	invalid_rho <- -1 / (d - 1)
	cor_model_invalid <- new("CSCorrelationModel", parameters = atanh(invalid_rho))
	expect_error(
		new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model_invalid),
		"Correlation rho of .* is not valid"
	)
})

test_that("Homogeneous-CS parameter setting and metadata work", {
	d <- 4L
	var_model <- new("HomogeneousVarianceModel", parameters = 0)
	cor_model <- new("CSCorrelationModel", parameters = 0)
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	new_params <- c(log(10), atanh(0.25))
	obj <- set_parameters(obj, new_params)
	expect_equal(get_parameters(obj), new_params)
	expect_error(set_parameters(obj, c(1, 2, 3)))

	expect_equal(n_parameters(obj), 2)
	expect_equal(get_start_values(obj), c(0, 0.1))
	expect_false(is_diagonal(obj))
	expect_equal(get_lower_bounds(obj), rep(-Inf, 2))

	obj_diag <- set_parameters(obj, c(1, 0))
	expect_true(is_diagonal(obj_diag))
})

test_that("Homogeneous-CS computations are correct", {
	d <- 3L
	st_dev <- 2
	rho <- 0.5
	var_model <- new("HomogeneousVarianceModel", parameters = log(st_dev^2))
	cor_model <- new("CSCorrelationModel", parameters = atanh(rho))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	interp_params <- get_interpretable_parameters(obj)
	expect_equal(interp_params$st_dev, st_dev)
	expect_equal(interp_params$correlation, rho)

	R <- matrix(rho, d, d); diag(R) <- 1
	Sigma_expected <- st_dev^2 * R

	expect_equal(as.matrix(compute_covariance_matrix(obj)), Sigma_expected)
	expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), as.matrix(solve(Sigma_expected)), tolerance = 1e-7)
	expect_equal(compute_log_det_covariance_matrix(obj), as.numeric(determinant(Sigma_expected, logarithm = TRUE)$modulus))
})

test_that("show method for Homogeneous-CS model works", {
	d <- 2L
	st_dev <- 2
	rho <- 0.5
	var_model <- new("HomogeneousVarianceModel", parameters = log(st_dev^2))
	cor_model <- new("CSCorrelationModel", parameters = atanh(rho))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	expect_output(show(obj), "Homogeneous-CS Covariance object \\(dimension: 2\\)")
	expect_output(show(obj), "  Standard Deviation \\(sigma\\): 2\\.0000")
	expect_output(show(obj), "  Correlation \\(rho\\): 0\\.5000")
})

# Tests for FlexibleCovariance with Heterogeneous CS Structure
# (Replaces the old HeterogeneousCSCov tests)

context("Testing FlexibleCovariance with Heterogeneous CS Structure")

test_that("Heterogeneous-CS computations are correct", {
	d <- 3L
	st_devs <- c(2, 3, 4)
	rho <- 0.5
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = log(st_devs^2))
	cor_model <- new("CSCorrelationModel", parameters = atanh(rho))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	interp_params <- get_interpretable_parameters(obj)
	expect_equal(interp_params$st_devs, st_devs)
	expect_equal(interp_params$correlation, rho)

	D <- diag(st_devs)
	R <- matrix(rho, d, d); diag(R) <- 1
	Sigma_expected <- D %*% R %*% D

	expect_equal(as.matrix(compute_covariance_matrix(obj)), Sigma_expected)
	expect_equal(as.matrix(get_cholesky_factor(obj)), chol(Sigma_expected))
	expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), solve(Sigma_expected))
	expect_equal(compute_log_det_covariance_matrix(obj), as.numeric(determinant(Sigma_expected)$modulus))
})

# Tests for FlexibleCovariance with AR1 Structure

context("Testing FlexibleCovariance with AR1 Structure")

test_that("Homogeneous-AR1 model computations are correct", {
	d <- 3L
	st_dev <- 2
	rho <- 0.5
	var_model <- new("HomogeneousVarianceModel", parameters = log(st_dev^2))
	cor_model <- new("AR1CorrelationModel", parameters = atanh(rho))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	time_diffs <- abs(outer(1:d, 1:d, "-"))
	R <- rho^time_diffs
	Sigma_expected <- st_dev^2 * R

	expect_equal(as.matrix(compute_covariance_matrix(obj)), Sigma_expected)
})

test_that("Heterogeneous-AR1 model computations are correct", {
	d <- 3L
	st_devs <- c(2, 3, 4)
	rho <- 0.5
	var_model <- new("HeterogeneousVarianceModel", dimension = d, parameters = log(st_devs^2))
	cor_model <- new("AR1CorrelationModel", parameters = atanh(rho))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	time_diffs <- abs(outer(1:d, 1:d, "-"))
	R <- rho^time_diffs
	D <- diag(st_devs)
	Sigma_expected <- D %*% R %*% D

	expect_equal(as.matrix(compute_covariance_matrix(obj)), Sigma_expected)
})

# Tests for FlexibleCovariance with Kernel Structures

context("Testing FlexibleCovariance with Kernel Structures")

test_that("Squared Exponential (Gaussian) kernel model works", {
	d <- 3L
	coords <- 1:d
	data_context <- list(coords = coords)
	length_scale <- 1.5
	st_dev <- 2

	var_model <- new("HomogeneousVarianceModel", parameters = log(st_dev^2))
	cor_model <- new("SquaredExpKernelModel", parameters = log(length_scale))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	dist_sq <- as.matrix(dist(coords)^2)
	R_expected <- exp(-dist_sq / (2 * length_scale^2))
	Sigma_expected <- st_dev^2 * R_expected

	expect_equal(as.matrix(compute_covariance_matrix(obj, data_context)), Sigma_expected)
	expect_output(show(obj), "Homogeneous-SquaredExp Covariance object \\(dimension: 3\\)")
	expect_output(show(obj), "  length_scale: 1\\.5000")
})

test_that("Matern kernel model works", {
	d <- 3L
	coords <- 1:d
	data_context <- list(coords = coords)
	length_scale <- 1.5
	nu <- 0.5 # For nu=0.5, Matern is equivalent to Exponential
	st_dev <- 2

	var_model <- new("HomogeneousVarianceModel", parameters = log(st_dev^2))
	cor_model <- new("MaternKernelModel", parameters = c(log(length_scale), log(nu)))
	obj <- new("FlexibleCovariance", dimension = d, variance_model = var_model, correlation_model = cor_model)

	dists <- as.matrix(dist(coords))
	R_expected <- exp(-dists / length_scale) # Exponential formula:
	Sigma_expected <- st_dev^2 * R_expected

	expect_equal(as.matrix(compute_covariance_matrix(obj, data_context)), Sigma_expected)
	expect_output(show(obj), "Homogeneous-Matern Covariance object \\(dimension: 3\\)")
	expect_output(show(obj), "  length_scale: 1\\.5000")
	expect_output(show(obj), "  nu: 0\\.5000")
})
