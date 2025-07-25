#
# Helper Functions for Covariance Structures


##' Forward CS transformation: theta -> rho 
cs_theta_to_rho <- function(theta,n) {
    a <- 1/(n-1)
    rho <- plogis(theta) * (1 + a) - a  
    return(rho)
} 

##' Inverse CS transformation: rho -> theta 
cs_rho_to_theta <- function(rho, n) {
    a <- 1 / (n - 1)
    theta <- qlogis((rho + a) / (1 + a))

    return(theta)
}

##' Forward AR1 transformation: theta -> rho
ar1_theta_to_rho <- function(theta) {
    rho <- theta / sqrt(1 + theta^2)
    
    return(rho)
}

##' Inverse AR1 transformation: rho -> theta
ar1_rho_to_theta <- function(rho) {
    theta <- rho / sqrt(1 - rho^2)

    return(theta)
}

##' @title Get Vech Diagonal Indices
##'
##' @description Internal helper to find the 1-based indices of the diagonal elements
##' within a 'vech' vector (the vectorized lower triangle of a matrix).
##'
##' @param d An integer, the dimension of the square matrix.
##' @return A numeric vector of integer indices for the diagonal elements.
##' @keywords internal
vech_diag_indices <- function(d) {
    if (d == 0) return(integer(0))
    if (d == 1) return(1L)
    1L + c(0, cumsum(d:2))
}


##' @title Get Vech Indices
##'
##' @description Internal helper to get row and column indices for a 'vech'
##' (column-major, lower-triangular) ordered vector.
##'
##' @param d An integer, the dimension of the square matrix.
##' @return A list with components 'i' (row indices) and 'j' (column indices).
##' @keywords internal
get_vech_indices <- function(d) {
	if (d == 1) return(list(i = 1L, j = 1L))
	list(
		i = unlist(lapply(1:d, function(j) j:d)),
		j = unlist(lapply(1:d, function(j) rep(j, d - j + 1)))
	)
}

##' @title Reconstruct Cholesky Factor from Parameters
##'
##' @description Internal helper function that reconstructs the dense
##' lower-triangular Cholesky factor `L` from an unconstrained parameter vector.
##'
##' @param param_vec The vector of unconstrained parameters.
##' @param d The dimension of the matrix.
##' @return A `dtrMatrix` object representing the Cholesky factor `L`.
##' @keywords internal
get_chol_from_params <- function(param_vec, d) {

	L_vec <- param_vec
	diag_indices <- vech_diag_indices(d)
	L_vec[diag_indices] <- exp(L_vec[diag_indices])

	L <- matrix(0, nrow = d, ncol = d)
	vech_indices <- get_vech_indices(d)
	L[cbind(vech_indices$i, vech_indices$j)] <- L_vec

	as(as(as(L, "dMatrix"), "triangularMatrix"), "unpackedMatrix")
}

##' @title Get Lower-Triangular Indices
##' @description Internal helper function computes the 1-based row and column
##'   indices for elements in the lower triangle (including the diagonal) of a
##'   square matrix. The indices are returned in column-major order.
##' @param d An integer, the dimension of the square matrix.
##' @return A list with two integer vector components: `i` (row indices)
##'   and `j` (column indices).
##' @keywords internal
get_lower_tri_indices <- function(d) {
    list(
        i = unlist(lapply(1:d, function(j) j:d)),
        j = unlist(lapply(1:d, function(j) rep(j, d - j + 1)))
    )
}

##' Force Matrix to Symmetric Dense Format
##' 
##' Internal helper to convert matrices to dsyMatrix class using modern Matrix package syntax.
##'
##' @param x A matrix-like object to convert
##' @return A dsyMatrix object
##' @keywords internal
force_dsyMatrix <- function(x) {
    as(as(as(x, "dMatrix"), "symmetricMatrix"), "unpackedMatrix")
}

##' Helper function for log determinant computation of structured covariance matrices
##'
##' Computes log determinant using decomposition: log(det(D*R*D)) = log(det(D²)) + log(det(R))
##' where D is variance scaling and R is correlation matrix.
##'
##' @param object A structured covariance object
##' @return Numeric log determinant value
##' @keywords internal
compute_log_det_structured <- function(object) {
    d <- object@dimension
      
    R <- compute_correlation_matrix(object)
    log_det_R <- as.numeric(Matrix::determinant(R, logarithm = TRUE)$modulus)
    
    log_det_V <- if (is(object, "HomogeneousVariance")) {
        # All variances equal: log(det(sigma^2 * I)) = d * log(sigma^2)
        d * object@vparameters[1]
    } else {
        # Individual variances: log(det(D^2)) = sum(log(variances))
        sum(object@vparameters)
    }
    
    return(log_det_V + log_det_R)
}

##' Helper function for inverse computation of structured covariance matrices
##'
##' Computes inverse using decomposition: inv(D*R*D) = inv(D)*inv(R)*inv(D)
##' where D is variance scaling and R is correlation matrix.
##'
##' @param object A structured covariance object
##' @return Inverse covariance matrix as dsyMatrix
##' @keywords internal
compute_inverse_structured <- function(object) {
    d <- object@dimension
    
    R <- compute_correlation_matrix(object)
    inv_R <- solve(R)
    
    if (is(object, "HomogeneousVariance")) {
        # inv(sigma^2 * R) = (1/sigma^2) * inv(R)
        inv_sigma_sq <- exp(-object@vparameters[1])
        inv_Sigma <- inv_sigma_sq * inv_R
    } else {
        # inv(D * R * D) = inv(D) * inv(R) * inv(D)
        inv_st_devs <- exp(-0.5 * object@vparameters)
        inv_D <- Diagonal(d, x = inv_st_devs)
        inv_Sigma <- inv_D %*% inv_R %*% inv_D
    }
    
    force_dsyMatrix(inv_Sigma)
}

# Helper Functions for CS and AR1 Covariance Structures

##' Get Start Values for Structured Covariance Models
##' 
##' @description Helper function to generate start values for CS and AR1 covariance structures.
##' @param object A covariance structure object (CS or AR1).
##' @return Numeric vector of start values: variance parameters (1.0) followed by 
##'   correlation parameter (0.0) if dimension > 1.
##' @keywords internal
get_structured_start_values <- function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    v_starts <- rep(0.0, n_v_params)
    
    if (d > 1) {
        c_starts <- 0.0
        return(c(v_starts, c_starts))
    } else {
        return(v_starts)
    }
}

##' Get Lower Bounds for Structured Covariance Models
##' 
##' @description Helper function to generate parameter bounds for CS and AR1 covariance structures.
##' @param object A covariance structure object (CS or AR1).
##' @return Numeric vector of lower bounds: variance bounds (0) followed by 
##'   correlation bound (-Inf) if dimension > 1.
##' @keywords internal
get_structured_lower_bounds <- function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    v_low <- rep(0, n_v_params)  # Variance bounds (log scale)
    
    if (d > 1) {
        c_low <- -Inf  # Correlation bound (atanh scale)
        return(c(v_low, c_low))
    } else {
        return(v_low)
    }
}

##' Create vech position to matrix distance mapping for structured covariance
##'
##' Maps vectorized half positions to matrix coordinates and distances for 
##' distance-based covariance structures like AR1.
##'
##' @param d Matrix dimension
##' @return List mapping vech positions to matrix coordinates, distances, and parameter indices
##' @keywords internal
get_vech_distance_mapping <- function(d) {
    # For dimension d, create mapping from vech position to matrix distance
    mapping <- list()
    vech_pos <- 1
    
    for (j in 1:d) {        # column
        for (i in j:d) {    # row (lower triangle)
            distance <- abs(i - j)
            mapping[[vech_pos]] <- list(
                matrix_pos = c(i, j),
                distance = distance,
                param_index = if (distance == 0) 1 else (1 + distance)  # 1=variance, 2+=correlations
            )
            vech_pos <- vech_pos + 1
        }
    }
    
    return(mapping)
}

##' Check if merMod Object Has Structured Covariance
##'
##' Utility function to determine if a fitted merMod object contains
##' structured covariance information.
##'
##' @param x A merMod object
##' @return Logical indicating presence of structure information
##' @export
has_structured_covariance <- function(x) {
    !is.null(attr(x, "cov_structures"))
}

##' Get Structure Types from merMod Object
##'
##' Utility function to extract the covariance structure types used in a model.
##'
##' @param x A merMod object  
##' @return Character vector of structure types, or NULL if unstructured
##' @export
get_structure_types <- function(x) {
    attr(x, "structure_types")
}
