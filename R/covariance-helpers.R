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
