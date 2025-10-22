
# Helper function to get a full suite of parameters for a given object
get_all_test_params <- function(obj) {
    d <- obj@dimension
    if (is(obj, "UnstructuredCovariance")) {
        n_v_params <- 0L  
    } else {
        n_v_params <- if (is(obj, "HomogeneousVariance")) 1L else d
    }

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
    expect_equal(default_params, rep(1, d))  
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
    expect_equal(get_start_values(obj), rep(1, d))
    expect_equal(get_lower_bounds(obj), rep(0, d))
    expect_equal(n_parameters(obj), d)

    # set & get garameters
    params <- c(2, 3, 4)
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
 
})

test_that("HeterogeneousDiagonalCovariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    
    theta <- c(2, 3, 4)
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    
    expected_sds <- theta
    expect_equal(lambdat_x_result, expected_sds)
     
    # edge case: dimension 1
    obj1 <- new("HeterogeneousDiagonalCovariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, 3)
    expect_equal(lambdat_x1, 3)   
})

test_that("HeterogeneousDiagonalCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("HeterogeneousDiagonalCovariance", dimension = d)
    obj <- set_parameters(obj, c(1, 2, 3))
    
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
    
    # small std devs
    small_params <- c(1e-3, 1e-3)
    obj_small <- set_parameters(obj, small_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))

    
    # large std dev
    large_params <- c(1e3, 1e3)
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
    expect_equal(default_params, 1)  

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
    
    expect_equal(get_start_values(obj), 1)  
    expect_equal(get_lower_bounds(obj), 0)  
 
    # set & get parameters
    param <- 3  
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
})

test_that("HomogeneousDiagonalCovariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
     
    theta <- 4
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    expected_sd <- theta 
    expect_equal(lambdat_x_result, expected_sd) 

    lind_result <- get_lind(obj)
    expect_equal(lind_result, rep(1L, d))
    
    # edge case: dimension 1
    obj1 <- new("HomogeneousDiagonalCovariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, 2)
    expect_equal(lambdat_x1, 2)
 })

test_that("HomogeneousDiagonalCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("HomogeneousDiagonalCovariance", dimension = d)
    obj <- set_parameters(obj, 2)
    
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
    
    # very std dev 
    small_param <- 1e-3
    obj_small <- set_parameters(obj, small_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large variance
    large_param <- 1e3
    obj_large <- set_parameters(obj, large_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
})


# UNSTRUCTURED COVARIANCE STRUCTURE TESTS

context("Unstructured Covariance Models Testing")

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
})

test_that("UnstructuredCovariance: method coverage", {
    d <- 2L
    obj <- new("UnstructuredCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    
    # object creation and validation
    expect_s4_class(obj, "UnstructuredCovariance")
    expect_true(validObject(obj))
    expected_params <- d * (d + 1) / 2  
    expect_equal(test_params$n_total, expected_params)
    expect_equal(test_params$n_v_params, 0)  
    expect_equal(test_params$n_c_params, expected_params)      
    
    ## start values
    expect_equal(get_start_values(obj), c(1, 0, 1))  

    # lower bounds 
    lower_bounds <- get_lower_bounds(obj)
    expect_equal(lower_bounds, c(0, -Inf, 0))
    
    # set & get parameters 
    params <- c(2, 0.5, sqrt(8.75))
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    
    # interpretable parameters
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_devs, c(2, 3))
    expect_equal(interpretable$correlation[1,2], 1/6)
    
    expect_true(is(get_lambda(obj), "sparseMatrix"))
    expect_equal(get_lind(obj), seq_len(expected_params))

    # covariance matrix
    expected_sigma <- matrix(c(4, 1, 1, 9), 2, 2)
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), expected_sigma)
    
    # cholesky factor
    computed_chol <- compute_cholesky_factor(obj)
    expect_equal(as.matrix(computed_chol), matrix(c(2, 0.5, 0, sqrt(8.75)), 2, 2))
    
    # inverse covariance matrix
    expected_inv <- solve(expected_sigma)
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))  
})

test_that("UnstructuredCovariance: compute_lambdat_x method", {
    d <- 2L
    obj <- new("UnstructuredCovariance", dimension = d)
    
    theta <- c(2, 0.5, 3)
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    
    expected_result <- c(2, 0.5, 3)
    expect_equal(lambdat_x_result, expected_result)
    
    # edge case: dimension 1
    obj1 <- new("UnstructuredCovariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, 3)
    expect_equal(lambdat_x1, 3)
})

test_that("UnstructuredCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("UnstructuredCovariance", dimension = d)
    params <- c(2, 0.3, 1.5, -0.2, 0.1, 1.8)
    obj <- set_parameters(obj, params)
    
    # computed matrices
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat))
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))
})

test_that("UnstructuredCovariance: boundary behavior", {
    d <- 2L
    obj <- new("UnstructuredCovariance", dimension = d)
    
    # very small diagonal elements
    small_params <- c(1e-3, 0.1, 1e-3)
    obj_small <- set_parameters(obj, small_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large diagonal elements
    large_params <- c(1e3, 0.1, 1e3)
    obj_large <- set_parameters(obj, large_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
    
    # near-singular case (small diagonal element)
    singular_params <- c(2, 1.9, 1e-6)  # L21 close to L11, L22 very small
    obj_singular <- set_parameters(obj, singular_params)
    cov_singular <- compute_covariance_matrix(obj_singular)
    
    # Should have very small determinant but still be computable
    expect_true(det(as.matrix(cov_singular)) < 1e-10)
})

# COMPOUND SYMMETRY (CS) COVARIANCE STRUCTURE TESTS

# HETEROGENEOUS CS COVARIANCE TESTS
test_that("HeterogeneousCSCovariance: object initialization", {
    d <- 3L
    obj <- new("HeterogeneousCSCovariance", dimension = d)
    
    expect_s4_class(obj, "HeterogeneousCSCovariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    expected_params <- d + 1  
    default_params <- get_parameters(obj)
    expect_equal(length(default_params), expected_params)
    expect_equal(default_params, get_start_values(obj)) 
})

test_that("HeterogeneousCSCovariance: method coverage", {
    d <- 3L
    obj <- new("HeterogeneousCSCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    
    # object creation and validation
    expect_s4_class(obj, "HeterogeneousCSCovariance")
    expect_true(validObject(obj))
    expected_params <- d + 1  
    expect_equal(test_params$n_total, expected_params)
    expect_equal(test_params$n_v_params, d)  
    expect_equal(test_params$n_c_params, 1)      
    
    # start values
    expect_equal(get_start_values(obj), c(rep(1.0, d), 0.0))  

    # lower bounds 
    lower_bounds <- get_lower_bounds(obj)
    expect_equal(lower_bounds, c(rep(0, d), -Inf))
    
    # set & get parameters 
    params <- c(2, 1.5, 3, 0.3)  
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    
    # interpretable parameters
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_devs, c(2, 1.5, 3))
    
    a <- 1/(d-1)  # a = 0.5 for d=3
    expected_rho <- plogis(0.3) * (1 + a) - a
    expect_equal(interpretable$correlation, expected_rho)

    lind_result <- get_lind(obj)
    expected_lind <- c(1L, 4L, 4L, 2L, 4L, 3L)
    expect_equal(lind_result,expected_lind)  

    # covariance matrix 
    R <- matrix(expected_rho, d, d)
    diag(R) <- 1.0
    D <- diag(c(2, 1.5, 3))
    expected_sigma <- D %*% R %*% D
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), as.matrix(expected_sigma))
    
    # cholesky factor
    computed_chol <- compute_cholesky_factor(obj)
    expected_chol <- t(chol(expected_sigma))
    expect_equal(as.matrix(computed_chol), as.matrix(expected_chol))
    
    # inverse covariance matrix
    expected_inv <- solve(expected_sigma)
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))

    # correlation matrix
    expected_cor <- compute_correlation_matrix(obj)
    computed_cor <- cov2cor(compute_covariance_matrix(obj))
    expect_equal(as.matrix(expected_cor), as.matrix(computed_cor))
})

test_that("HeterogeneousCSCovariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HeterogeneousCSCovariance", dimension = d)
    
    theta <- c(2, 1.5, 3, 0.2)  
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    
    expected_length <- d + 1
    expect_equal(length(lambdat_x_result), expected_length)
    
    # verify against set_parameters + compute_cholesky_factor
    obj_set <- set_parameters(obj, theta)
    L <- compute_cholesky_factor(obj_set)
    expected_from_full <- L[lower.tri(L, diag = TRUE)]
    skip()
    expect_equal(lambdat_x_result, as.numeric(expected_from_full))
})

test_that("HeterogeneousCSCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("HeterogeneousCSCovariance", dimension = d)
    params <- c(1, 2, 3, 0.4)  
    obj <- set_parameters(obj, params)
    
    # computed matrices
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat))
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))
    
})

test_that("HeterogeneousCSCovariance: boundary behavior", {
    d <- 2L
    obj <- new("HeterogeneousCSCovariance", dimension = d)
    
    # very small standard deviations
    small_params <- c(1e-3, 1e-3, 0.1)
    obj_small <- set_parameters(obj, small_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large standard deviations
    large_params <- c(1e3, 1e3, -0.2)
    obj_large <- set_parameters(obj, large_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
    
})

# HOMOGENEOUS CS COVARIANCE TESTS

test_that("HomogeneousCSCovariance: object initialization", {
    d <- 4L
    obj <- new("HomogeneousCSCovariance", dimension = d)
    
    expect_s4_class(obj, "HomogeneousCSCovariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    expected_params <- if (d > 1) 2L else 1L  
    default_params <- get_parameters(obj)
    expect_equal(length(default_params), expected_params)
    expect_equal(default_params, get_start_values(obj))
})

test_that("HomogeneousCSCovariance: method coverage", {
    d <- 4L
    obj <- new("HomogeneousCSCovariance", dimension = d)
    test_params <- get_all_test_params(obj)
    
    # object creation and validation
    expect_s4_class(obj, "HomogeneousCSCovariance")
    expect_true(validObject(obj))
    expect_equal(test_params$n_total, 2)
    expect_equal(test_params$n_v_params, 1)
    expect_equal(test_params$n_c_params, 1)
    
    expect_equal(get_start_values(obj), c(1.0, 0.0))
    expect_equal(get_lower_bounds(obj), c(0, -Inf))
    
    # set & get parameters 
    params <- c(2, 0)  
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_dev, 2)
    
    # CS correlation
    a <- 1/(d-1)
    expected_rho <- plogis(0) * (1 + a) - a
    expect_equal(interpretable$correlation, expected_rho)
    
    lind_result <- get_lind(obj)
    expected_lind <- c(1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L)
    expect_equal(lind_result, expected_lind)

    # covariance matrix 
    R <- matrix(expected_rho, d, d)
    diag(R) <- 1.0
    expected_sigma <- 2^2 * R  
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), as.matrix(expected_sigma))
    
    # cholesky factor
    expected_chol <- t(chol(expected_sigma))
    computed_chol <- compute_cholesky_factor(obj)
    expect_equal(as.matrix(computed_chol), as.matrix(expected_chol))
    
    # inverse covariance matrix
    expected_inv <- solve(expected_sigma)
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))

    # correlation matrix
    expected_cor <- compute_correlation_matrix(obj)
    computed_cor <- cov2cor(compute_covariance_matrix(obj))
    expect_equal(as.matrix(expected_cor), as.matrix(computed_cor))

})

test_that("HomogeneousCSCovariance: compute_lambdat_x method", { 
    d <- 3L  
    obj <- new("HomogeneousCSCovariance", dimension = d)
    theta <- c(2, 0)  
    lambdat_x_result <- compute_lambdat_x(obj, theta)

    a <- 1/(d-1)
    rho <- plogis(0) * (1 + a) - a
    R <- matrix(rho, d, d)
    diag(R) <- 1.0
    expected_Sigma <- 4 * R  
    expected_L <- t(chol(expected_Sigma))
    expected_vech <- expected_L[lower.tri(expected_L, diag = TRUE)]

    skip()
    expect_equal(lambdat_x_result, as.numeric(expected_vech))
        
    # edge case: dimension 1
    obj1 <- new("HomogeneousCSCovariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, 2)
    expect_equal(lambdat_x1, 2)
})

test_that("HomogeneousCSCovariance: mathematical consistency", {
    d <- 3L
    obj <- new("HomogeneousCSCovariance", dimension = d)
    obj <- set_parameters(obj, c(2, 0.3))  

    # computed matrices
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat))
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))    
})

test_that("HomogeneousCSCovariance: boundary behavior", {
    d <- 2L
    obj <- new("HomogeneousCSCovariance", dimension = d)
    
    # very small std dev
    small_param <- c(1e-3, 0.1)
    obj_small <- set_parameters(obj, small_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large std dev
    large_param <- c(1e3, -0.2)
    obj_large <- set_parameters(obj, large_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
})

# AR1 COVARIANCE STRUCTURE TESTS

# HETEROGENEOUS AR1 COVARIANCE TESTS

test_that("HeterogeneousAR1Covariance: object initialization", {
    d <- 3L
    obj <- new("HeterogeneousAR1Covariance", dimension = d)
    
    expect_s4_class(obj, "HeterogeneousAR1Covariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    expected_params <- d + 1 
    default_params <- get_parameters(obj)
    expect_equal(length(default_params), expected_params)
    expect_equal(default_params, get_start_values(obj))
    
})

test_that("HeterogeneousAR1Covariance: method coverage", {
    d <- 3L
    obj <- new("HeterogeneousAR1Covariance", dimension = d)
    test_params <- get_all_test_params(obj)
    
    # object creation and validation
    expect_s4_class(obj, "HeterogeneousAR1Covariance")
    expect_true(validObject(obj))
    expected_params <- d + 1  
    expect_equal(test_params$n_total, expected_params)
    expect_equal(test_params$n_v_params, d)  
    expect_equal(test_params$n_c_params, 1)      
    
    # start values
    expect_equal(get_start_values(obj), c(rep(1.0, d), 0.0))  

    # lower bounds 
    lower_bounds <- get_lower_bounds(obj)
    expect_equal(lower_bounds, c(rep(0, d), -Inf))
    
    # set & get parameters 
    params <- c(2, 1.5, 3, 0.5)  
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    
    # interpretable parameters
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_devs, c(2, 1.5, 3))
    
    # AR1 correlation transformation
    expected_rho <- ar1_theta_to_rho(0.5)  
    expect_equal(interpretable$correlation, expected_rho)
    
    
    # test get_lind 
    lind_result <- get_lind(obj)
    expected_length <- d * (d + 1) / 2
    expected_lind <- c(1 ,4 ,4 ,2 ,4 ,3)
    expect_equal(lind_result, expected_lind)  

    # covariance matrix 
    rho <- expected_rho
    R <- matrix(0, d, d)
    for (i in 1:d) {
        for (j in 1:d) {
            R[i,j] <- rho^abs(i-j)
        }
    }
    D <- diag(c(2, 1.5, 3))
    expected_sigma <- D %*% R %*% D
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), as.matrix(expected_sigma))
    
    # cholesky factor
    computed_chol <- compute_cholesky_factor(obj)
    expected_chol <- t(chol(expected_sigma))
    expect_equal(as.matrix(computed_chol), as.matrix(expected_chol))
    
    # inverse covariance matrix
    expected_inv <- solve(expected_sigma)
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))

    # correlation matrix
    expected_cor <- compute_correlation_matrix(obj)
    computed_cor <- cov2cor(compute_covariance_matrix(obj))
    expect_equal(as.matrix(expected_cor), as.matrix(computed_cor))
})

test_that("HeterogeneousAR1Covariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HeterogeneousAR1Covariance", dimension = d)
    
    theta <- c(2, 1.5, 3, 0.3)  
    lambdat_x_result <- compute_lambdat_x(obj, theta)
    
    expected_length <- d * (d + 1) / 2
    expect_equal(length(lambdat_x_result), expected_length)
    
    # verify against set_parameters + compute_cholesky_factor
    obj_set <- set_parameters(obj, theta)
    L <- compute_cholesky_factor(obj_set)
    expected_from_full <- L[lower.tri(L, diag = TRUE)]
    expect_equal(lambdat_x_result, as.numeric(expected_from_full))
    
    # edge case: dimension 1
    obj1 <- new("HeterogeneousAR1Covariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, 3)
    skip()
    expect_equal(lambdat_x1, 3)
})

test_that("HeterogeneousAR1Covariance: mathematical consistency", {
    d <- 3L
    obj <- new("HeterogeneousAR1Covariance", dimension = d)
    params <- c(1, 2, 3, 0.4)
    obj <- set_parameters(obj, params)
    
    # computed matrices
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat))
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))
})

test_that("HeterogeneousAR1Covariance: boundary behavior", {
    d <- 2L
    obj <- new("HeterogeneousAR1Covariance", dimension = d)
    
    # very small standard deviations
    small_params <- c(1e-3, 1e-3, 0.1)
    obj_small <- set_parameters(obj, small_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large standard deviations
    large_params <- c(1e3, 1e3, -0.2)
    obj_large <- set_parameters(obj, large_params)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
    
    # extreme correlation parameters
    extreme_params <- c(1, 1, -5)  
    obj_extreme <- set_parameters(obj, extreme_params)
    rho_extreme <- get_interpretable_parameters(obj_extreme)$correlation
    expect_true(rho_extreme >= -1.0 && rho_extreme <= 1.0)
})

# HOMOGENEOUS AR1 COVARIANCE TESTS

test_that("HomogeneousAR1Covariance: object initialization", {
    d <- 4L
    obj <- new("HomogeneousAR1Covariance", dimension = d)
    
    expect_s4_class(obj, "HomogeneousAR1Covariance")
    expect_true(validObject(obj))
    expect_equal(obj@dimension, d)
    
    expected_params <- if (d > 1) 2L else 1L  
    default_params <- get_parameters(obj)
    expect_equal(length(default_params), expected_params)
    expect_equal(default_params, get_start_values(obj))
    
})

test_that("HomogeneousAR1Covariance: method coverage", {
    d <- 4L
    obj <- new("HomogeneousAR1Covariance", dimension = d)
    test_params <- get_all_test_params(obj)
    
    # object creation and validation
    expect_s4_class(obj, "HomogeneousAR1Covariance")
    expect_true(validObject(obj))
    expect_equal(test_params$n_total, 2)
    expect_equal(test_params$n_v_params, 1)
    expect_equal(test_params$n_c_params, 1)
    
    expect_equal(get_start_values(obj), c(1.0, 0.0))
    expect_equal(get_lower_bounds(obj), c(0, -Inf))
    
    # set & get parameters 
    params <- c(2, 0)  
    obj <- set_parameters(obj, params)
    expect_equal(get_parameters(obj), params)
    interpretable <- get_interpretable_parameters(obj)
    expect_equal(interpretable$st_dev, 2)
    
    # correlation
    expected_rho <- ar1_theta_to_rho(0)  
    expect_equal(interpretable$correlation, expected_rho)
        
    # get_lind 
    lind_result <- get_lind(obj)
    expected_length <- d * (d + 1) / 2
    expect_equal(length(lind_result), expected_length)
    # expect_equal(lind_result, c(1, 2, 2, 1, 2, 1))
    
    # covariance matrix 
    rho <- expected_rho  
    R <- matrix(0, d, d)
    for (i in 1:d) {
        for (j in 1:d) {
            R[i,j] <- rho^abs(i-j)
        }
    }
    expected_sigma <- 2^2 * R  
    computed_sigma <- compute_covariance_matrix(obj)
    expect_equal(as.matrix(computed_sigma), as.matrix(expected_sigma))
    
    # cholesky factor
    expected_chol <- t(chol(expected_sigma))
    computed_chol <- compute_cholesky_factor(obj)
    expect_equal(as.matrix(computed_chol), as.matrix(expected_chol))
    
    # inverse covariance matrix
    expected_inv <- solve(expected_sigma)
    computed_inv <- compute_inverse_covariance_matrix(obj)
    expect_equal(as.matrix(computed_inv), as.matrix(expected_inv))
    
    # correlation matrix
    expected_cor <- compute_correlation_matrix(obj)
    computed_cor <- cov2cor(compute_covariance_matrix(obj))
    expect_equal(as.matrix(expected_cor), as.matrix(computed_cor))
})

test_that("HomogeneousAR1Covariance: compute_lambdat_x method", {
    d <- 3L
    obj <- new("HomogeneousAR1Covariance", dimension = d)
    
    theta <- c(2, 0) 
    lambdat_x_result <- compute_lambdat_x(obj, theta)
        
    # verify against set_parameters + compute_cholesky_factor
    obj_set <- set_parameters(obj, theta)
    L <- compute_cholesky_factor(obj_set)
    expected_from_full <- L[lower.tri(L, diag = TRUE)]
    expect_equal(lambdat_x_result, as.numeric(expected_from_full))
    
    # edge case: dimension 1
    obj1 <- new("HomogeneousAR1Covariance", dimension = 1L)
    lambdat_x1 <- compute_lambdat_x(obj1, 2)
    expect_equal(lambdat_x1, 2)
})

test_that("HomogeneousAR1Covariance: mathematical consistency", {
    d <- 3L
    obj <- new("HomogeneousAR1Covariance", dimension = d)
    obj <- set_parameters(obj, c(2, 0.3))  

    # computed matrices
    cov_mat <- compute_covariance_matrix(obj)
    chol_mat <- compute_cholesky_factor(obj)
    inv_mat <- compute_inverse_covariance_matrix(obj)
    
    # L * L' = Sigma
    reconstructed_sigma <- chol_mat %*% t(chol_mat)
    expect_equal(as.matrix(reconstructed_sigma), as.matrix(cov_mat))
    
    # Sigma * Sigma^(-1) = I
    identity_check <- cov_mat %*% inv_mat
    expect_equal(as.matrix(identity_check), diag(d))
})

test_that("HomogeneousAR1Covariance: boundary behavior", {
    d <- 2L
    obj <- new("HomogeneousAR1Covariance", dimension = d)
    
    # very small std dev
    small_param <- c(1e-3, 0.1)
    obj_small <- set_parameters(obj, small_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_small)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_small)))))
    
    # large std dev
    large_param <- c(1e3, -0.2)
    obj_large <- set_parameters(obj, large_param)
    
    expect_true(all(is.finite(as.matrix(compute_covariance_matrix(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_cholesky_factor(obj_large)))))
    expect_true(all(is.finite(as.matrix(compute_inverse_covariance_matrix(obj_large)))))
})


