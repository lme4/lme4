#
# Helper Functions for Covariance Structures (Refactored)


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
	if (d == 0) return(list(i = integer(0), j = integer(0)))
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
	if (d == 0) return(new("dtrMatrix", uplo = "L", Dim = c(0L, 0L)))

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
    if (d == 0) return(list(i = integer(0), j = integer(0)))
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
    if (d == 0) return(0)
    
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
    if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
    
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
