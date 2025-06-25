library(testthat)
library(Matrix)

# Load the S4 class definitions and helper functions.
# source("~/Documents/R/lme4/R/covariance_structures_refactored.R")
# source("~/Documents/R/lme4/R/covariance-helpers.R")

# Helper function to get a full suite of parameters for a given object
get_all_test_params <- function(obj) {
    d <- obj@dimension
    n_v_params <- if (is(obj, "HomogeneousVariance")) 1L else d

    list(
        n_v_params = n_v_params,
        n_c_params = n_parameters(obj) - n_v_params,
        n_total = n_parameters(obj)
    )
}

# Testing Diagonal Covariance Models 

context("Testing Diagonal Covariance Models")

test_that("HeterogeneousDiagonalCovariance: all methods work", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    test_params <- get_all_test_params(obj)

    # Creation & Parameters
    expect_s4_class(obj, "HeterogeneousDiagonalCovariance")
    expect_true(validObject(obj))
    expect_equal(test_params$n_total, d)
    expect_equal(get_start_values(obj), rep(0, d))

    # Set & Get
    params <- log(c(2^2, 3^2, 4^2))
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    expect_equal(get_interpretable_parameters(obj)$st_devs, c(2, 3, 4))
    expect_true(is_diagonal(obj))
    expect_equal(get_lower_bounds(obj), rep(-Inf, d))

    # Computation
    expect_equal(as.matrix(get_lambda(obj)), as.matrix(Diagonal(d)))
    expect_equal(get_lind(obj), integer(0))

    Sigma <- Diagonal(d, x = c(4, 9, 16))
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))
    expect_equal(as.matrix(get_cholesky_factor(obj)), as.matrix(Diagonal(d, x = c(2, 3, 4))))
    expect_equal(compute_log_det_covariance_matrix(obj), sum(params))
    expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), as.matrix(Diagonal(d, x = 1 / c(4, 9, 16))))

    # Show 
    expect_output(show(obj), "HeterogeneousDiagonalCovariance object \\(dimension: 3, parameters: 3\\)")
    expect_output(show(obj), "Standard Deviations: 2.0000 3.0000 4.0000")
})

test_that("HomogeneousDiagonalCovariance: all methods work", {
    d <- 4L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
    test_params <- get_all_test_params(obj)

    # Creation & Parameters
    expect_s4_class(obj, "HomogeneousDiagonalCovariance")
    expect_equal(test_params$n_total, 1)
    expect_equal(get_lower_bounds(obj), -Inf)

    # Set & Get
    obj <- set_parameters(obj, log(3^2))
    expect_equal(get_interpretable_parameters(obj)$st_dev, 3)

    # Computation
    Sigma <- Diagonal(d, x = rep(9, d))
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))
    expect_equal(compute_log_det_covariance_matrix(obj), 4 * log(3^2))
    
    # Show
    output_lines <- capture.output(show(obj))
    expect_match(output_lines[1], "HomogeneousDiagonalCovariance object \\(dimension: 4, parameters: 1\\)")
    expect_match(output_lines[2], "Interpretable Parameters:")
    expect_match(output_lines[3], "Standard Deviation \\(sigma\\): 3.0000")
})


# Testing Unstructured Covariance Models

context("Testing Unstructured Covariance Models")

test_that("HeterogeneousUnstructuredCovariance: all methods work", {
    d <- 2L
    obj <- new("HeterogeneousUnstructuredCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    expect_equal(test_params$n_total, d + d * (d + 1) / 2) # 2 + 3 = 5

    # Set & Get
    v_params <- log(c(2^2, 3^2))
    c_params <- c(log(1.5), 0.5, log(2.5)) # L = [1.5, 0; 0.5, 2.5]
     obj <- set_parameters(obj, c(v_params, c_params))
    expect_false(is_diagonal(obj))
    L <- matrix(c(1.5, 0.5, 0, 2.5), 2, 2)
    Sigma <- tcrossprod(L)
    expect_equal(get_interpretable_parameters(obj)$st_devs, sqrt(diag(Sigma)))
    
    # Computation
    expect_true(is(get_lambda(obj), "sparseMatrix"))
    expect_equal(get_lind(obj), (test_params$n_v_params + 1):test_params$n_total)

    expect_equal(as.matrix(get_cholesky_factor(obj)), L)
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))
    expect_equal(compute_log_det_covariance_matrix(obj), 2 * sum(c(log(1.5), log(2.5))))
    expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), as.matrix(solve(Sigma)))
})

test_that("HomogeneousUnstructuredCovariance: all methods work", {
    d <- 2L
    obj <- new("HomogeneousUnstructuredCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    expect_equal(test_params$n_total, 1 + d * (d + 1) / 2) # 1 + 3 = 4

    # Set & Get
    v_param <- log(5^2)
    c_params <- c(log(1.5), 0.5, log(2.5))
    obj <- set_parameters(obj, c(v_param, c_params))
    L <- matrix(c(1.5, 0.5, 0, 2.5), 2, 2)
    Sigma <- tcrossprod(L)
    expect_equal(get_interpretable_parameters(obj)$st_dev, mean(sqrt(diag(Sigma))))

    expect_equal(get_lind(obj), (test_params$n_v_params + 1):test_params$n_total)

    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))

    expect_output(show(obj), "HomogeneousUnstructuredCovariance object \\(dimension: 2, parameters: 4\\)")
})

# Testing Compound Symmetry (CS) Covariance Models

context("Testing CS Covariance Models")

test_that("HeterogeneousCSCovariance: all methods work", {
    d <- 3L
    obj <- new("HeterogeneousCSCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    expect_equal(test_params$n_total, d + 2) # 3 + 2 = 5

    # Set & Get
    v_params <- log(c(1.5^2, 2^2, 2.5^2))
    c_params <- c(0.1, atanh(0.6))
    obj <- set_parameters(obj, c(v_params, c_params))
    expect_equal(get_interpretable_parameters(obj)$correlation, 0.6)
    expect_equal(get_lower_bounds(obj), rep(-Inf, 5))

    # Computation 
    expect_equal(get_lind(obj), c(4, 5, 5, 4, 5, 4))
    expect_true(is(get_lambda(obj), "sparseMatrix"))

    D <- Diagonal(d, x = c(1.5, 2, 2.5))
    R <- matrix(0.6, d, d); diag(R) <- 1
    Sigma <- D %*% R %*% D
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))
    expect_equal(compute_log_det_covariance_matrix(obj), as.numeric(determinant(Sigma)$modulus))
    expect_equal(as.matrix(compute_inverse_covariance_matrix(obj)), as.matrix(solve(Sigma)))
    expect_equal(as.matrix(get_cholesky_factor(obj)), as.matrix(chol(Sigma)))

    # Show
    expect_output(show(obj), "HeterogeneousCSCovariance object \\(dimension: 3, parameters: 5\\)")
})

test_that("HomogeneousCSCovariance: all methods work", {
    d <- 2L
    obj <- new("HomogeneousCSCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    expect_equal(test_params$n_total, 1 + 2) # 1 + 2 = 3

    # Set & Get
    obj <- set_parameters(obj, c(log(4^2), 0.1, atanh(0.5)))
    expect_equal(get_interpretable_parameters(obj)$st_dev, 4)
    expect_equal(get_interpretable_parameters(obj)$correlation, 0.5)

    # Computation
    R <- matrix(0.5, d, d); diag(R) <- 1
    Sigma <- 4^2 * R
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))

    # Show
    expect_output(show(obj), "HomogeneousCSCovariance object \\(dimension: 2, parameters: 3\\)")
    expect_output(show(obj), "Standard Deviation \\(sigma\\): 4.0000")
    expect_output(show(obj), "Correlation \\(rho\\): 0.5000")
})

# Testing Autoregressive (AR1) Covariance Models 

context("Testing AR1 Covariance Models")

test_that("HeterogeneousAR1Covariance: all methods work", {
    d <- 3L
    obj <- new("HeterogeneousAR1Covariance", dimension = d)
    test_params <- get_all_test_params(obj)
    expect_equal(test_params$n_total, d + d) # 3 + 3 = 6

    # Set & Get
    v_params <- log(c(2^2, 3^2, 4^2))
    c_params <- c(0, atanh(0.7), 0.1) # lag0, lag1(rho), lag2
    obj <- set_parameters(obj, c(v_params, c_params))
    expect_equal(get_interpretable_parameters(obj)$st_devs, c(2, 3, 4))
    expect_equal(get_interpretable_parameters(obj)$correlation_lag1, 0.7)

    # Computation 
    expect_equal(get_lind(obj), c(4, 5, 6, 4, 5, 4))
    expect_true(is(get_lambda(obj), "sparseMatrix"))

    D <- Diagonal(d, x = c(2, 3, 4))
    R <- 0.7^abs(outer(1:d, 1:d, "-"))
    Sigma <- D %*% R %*% D
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))
    expect_equal(as.matrix(get_cholesky_factor(obj)), as.matrix(chol(Sigma)))

    # Show
    expect_output(show(obj), "HeterogeneousAR1Covariance object \\(dimension: 3, parameters: 6\\)")
})

test_that("HomogeneousAR1Covariance: all methods work", {
    d <- 4L
    obj <- new("HomogeneousAR1Covariance", dimension = d)
    test_params <- get_all_test_params(obj)
    expect_equal(test_params$n_total, 1 + d) # 1 + 4 = 5

    # Set & Get
    v_param <- log(3^2)
    c_params <- c(0, atanh(0.8), 0.1, 0.05)
    obj <- set_parameters(obj, c(v_param, c_params))
    expect_equal(get_interpretable_parameters(obj)$correlation_lag1, 0.8)

    # Computation
    indices <- get_lower_tri_indices(d)
    lag <- indices$i - indices$j
    expect_equal(get_lind(obj), test_params$n_v_params + lag + 1)
    expect_true(is(get_lambda(obj), "sparseMatrix"))

    
    R <- 0.8^abs(outer(1:d, 1:d, "-"))
    Sigma <- 3^2 * R
    expect_equal(as.matrix(compute_covariance_matrix(obj)), as.matrix(Sigma))
    expect_equal(compute_log_det_covariance_matrix(obj), as.numeric(determinant(Sigma)$modulus))
})


# --- Testing Edge Cases ---

context("Testing Edge Cases (Dimension 0 and 1)")

test_that("Dimension 0 objects are created correctly for all structures", {
    classes <- c("Diagonal", "Unstructured", "CS", "AR1")
    for (struct in classes) {
        class_name <- paste0("Homogeneous", struct, "Covariance")
        obj <- new(class_name, dimension = 0L)
        expect_true(validObject(obj))
        expect_equal(nrow(compute_covariance_matrix(obj)), 0)
        expect_equal(n_parameters(obj), 1)
  }
})

test_that("Dimension 1 objects are created correctly for all structures", {
    classes <- c("Diagonal", "Unstructured", "CS", "AR1")
    for (struct in classes) {
        class_name <- paste0("Heterogeneous", struct, "Covariance")
        obj <- new(class_name, dimension = 1L)
        expect_true(validObject(obj))
        obj <- set_parameters(obj, get_start_values(obj))
        expect_equal(as.matrix(compute_covariance_matrix(obj)), matrix(1))
        expect_true(Matrix::isDiagonal(build_correlation_matrix(obj)))
    }
})
