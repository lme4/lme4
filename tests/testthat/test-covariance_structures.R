
skip()

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

context("Diagonal Covariance Models Testing")

# HETEROGENEOUS DIAGONAL COVARIANCE 


test_that("HeterogeneousDiagonalCovariance: object initialization", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    
    expect_s4_class(obj, "HeterogeneousDiagonalCovariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    default_params <- get_parameters(obj)
    expect_equal(default_params, rep(0, d))  

    # test multiple dimensions
    for (test_d in c(1L, 2L, 5L, 10L)) {
        obj_multi <- new("HeterogeneousDiagonalCovariance", dimension = test_d)
        expect_true(validObject(obj_multi))
        expect_equal(n_parameters(obj_multi), test_d)
        expect_equal(length(get_parameters(obj_multi)), test_d)
    }
})

test_that("HeterogeneousDiagonalCovariance: method coverage", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    test_params <- get_all_test_params(obj)

    # object creation and validation 
    expect_s4_class(obj, "HeterogeneousDiagonalCovariance")
    expect_true(validObject(obj))
    expect_equal(test_params$n_total, d)
    expect_equal(test_params$n_v_params, d)
    expect_equal(test_params$n_c_params, 0) 
    expect_equal(get_start_values(obj), rep(0, d))
    expect_equal(get_lower_bounds(obj), rep(0, d))
    expect_equal(n_parameters(obj), d)

    # set & get garameters
    params <- log(c(4, 9, 16))  
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    expect_true(validObject(obj))
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_devs, c(2, 3, 4))
    expect_equal(as.matrix(get_lambda(obj)), as.matrix(Diagonal(d)))
    expect_equal(get_lind(obj), seq_len(d))


    # covariance matrix
    expected_sigma <- Diagonal(d, x = c(4, 9, 16))
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), as.matrix(expected_sigma))
    
    # cholesky factor
    expected_chol <- Diagonal(d, x = c(2, 3, 4))
    computed_chol <- compute_cholesky_factor(obj)
    expect_equal(as.matrix(computed_chol), as.matrix(expected_chol))
    
    # inverse covariance matrix
    expected_inv <- Diagonal(d, x = 1 / c(4, 9, 16))
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))

    # matrix type Consistency
    expect_true(is(computed_sigma, "Matrix"))  
    expect_true(is(computed_chol, "Matrix"))
    expect_true(is(computed_inv, "Matrix"))
 
})

test_that("HeterogeneousDiagonalCovariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    
    theta <- log(c(4, 9, 16))  # should this be log?
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    
    expected_variances <- exp(theta) 
    expect_equal(lambdat_x_result, expected_variance)
     
    # edge case: dimension 1
    obj1 <- new("HeterogeneousDiagonalCovariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, log(9))
    expect_equal(lambdat_x1, 9)   
})

test_that("HeterogeneousDiagonalCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    obj <- set_parameters(obj, log(c(1, 4, 9)))  
    
    # computed matrices
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    expect_equal(dim(cov_mat), c(d, d))
    expect_equal(dim(chol_mat), c(d, d))
    expect_equal(dim(inv_mat), c(d, d))
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat)) 
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))   
})

test_that("HeterogeneousDiagonalCovariance: boundary behavior", {
    d <- 2L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    
    # very small variances
    small_params <- log(c(1e-6, 1e-5))
    obj_small <- set_parameters(obj, small_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    expect_true(is.finite(compute_log_det_covariance_matrix(obj_small)))
    
    # large variances
    large_params <- log(c(1e6, 1e7))
    obj_large <- set_parameters(obj, large_params) 

    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
})

# HOMOGENEOUS DIAGONAL COVARIANCE 

test_that("HomogeneousDiagonalCovariance: object initialization", {
    d <- 4L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
    
    expect_s4_class(obj, "HomogeneousDiagonalCovariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    default_params <- get_parameters(obj)
    expect_equal(default_params, 0)  

    # test multiple dimensions
    for (test_d in c(1L, 2L, 5L, 10L)) {
        obj_multi <- new("HomogeneousDiagonalCovariance", dimension = test_d)
        expect_true(validObject(obj_multi))
        expect_equal(n_parameters(obj_multi), 1)  
    }
})

test_that("HomogeneousDiagonalCovariance: method coverage", {
    d <- 4L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    
    # object creation and validation 
    expect_s4_class(obj, "HomogeneousDiagonalCovariance")
    expect_true(validObject(obj))
    expect_equal(test_params$n_total, 1)
    expect_equal(test_params$n_v_params, 1)
    expect_equal(test_params$n_c_params, 0)  
    
    expect_equal(get_start_values(obj), 0)  
    expect_equal(get_lower_bounds(obj), 0)  
 
    # set & get parameters
    param <- log(9)  
    obj <- set_parameters(obj, param)
    expect_equal(get_parameters(obj), param)
    
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_dev, 3)
    expect_true(is_diagonal(obj))
 
    expect_equal(as.matrix(get_lambda(obj)), as.matrix(Diagonal(d)))
    expect_equal(get_lind(obj), rep(1L, d))  

    # covariance Matrix
    expected_sigma <- Diagonal(d, x = rep(9, d))  
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), as.matrix(expected_sigma))

    # cholesky Factor
    expected_chol <- Diagonal(d, x = rep(3, d)) 
    computed_chol <- compute_cholesky_factor(obj)
    expect_equal(as.matrix(computed_chol), as.matrix(expected_chol))
    
    # inverse Covariance Matrix
    expected_inv <- Diagonal(d, x = rep(1/9, d))  
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))
    
    # matrix type Consistency
    expect_true(is(computed_sigma, "Matrix"))  
    expect_true(is(computed_chol, "Matrix"))
    expect_true(is(computed_inv, "Matrix"))
})

test_that("HomogeneousDiagonalCovariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
     
    theta <- log(16)
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    expected_variance <- exp(theta)  
    expect_equal(lambdat_x_result, expected_variance) 

    lind_result <- get_lind(obj)
    expect_equal(lind_result, rep(1L, d))
    
    # edge case: dimension 1
    obj1 <- new("HomogeneousDiagonalCovariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, log(4))
    expect_equal(lambdat_x1, 4)
 })

test_that("HomogeneousDiagonalCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
    obj <- set_parameters(obj, log(2^2))  # variance = 4, std dev = 2
    
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    
    expect_equal(dim(cov_mat), c(d, d))
    expect_equal(dim(chol_mat), c(d, d))
    expect_equal(dim(inv_mat), c(d, d))
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat))
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))     
})

test_that("HomogeneousDiagonalCovariance: boundary behavior", {
    d <- 2L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
    
    # very small variance 
    small_param <- log(1e-6)
    obj_small <- set_parameters(obj, small_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large variance
    large_param <- log(1e6)
    obj_large <- set_parameters(obj, large_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
})


# UNSTRUCTURED COVARIANCE STRUCTURE TESTS

context("Comprehensive Unstructured Covariance Models Testing")

test_that("UnstructuredCovariance: object initialization", {
    d <- 3L
    obj <- new("UnstructuredCovariance", dimension = d)
    
    expect_s4_class(obj, "UnstructuredCovariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    expected_params <- d * (d + 1) / 2 
    default_params <- get_parameters(obj)
    expect_equal(length(default_params), expected_params)
    expect_equal(default_params, get_start_values(obj)) 
    
    # test multiple dimensions
    for (test_d in c(1L, 2L, 4L, 5L)) {
        obj_multi <- new("UnstructuredCovariance", dimension = test_d)
        expect_true(validObject(obj_multi))
        expected_multi_params <- test_d * (test_d + 1) / 2
        expect_equal(n_parameters(obj_multi), expected_multi_params)
        expect_equal(length(get_parameters(obj_multi)), expected_multi_params)
    }
})

 
