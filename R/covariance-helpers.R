##' @title Get Vech Diagonal Indices
##'
##' @description Internal helper to find the 1-based indices of the diagonal elements
##' within a 'vech' vector (the vectorized lower triangle of a matrix,
##' stacked column-wise).
##'
##' This is useful for selectively applying transformations (e.g., log/exp)
##' to the diagonal elements of a Cholesky factor's parameter vector.
##'
##' @param d An integer, the dimension of the square matrix.
##' @return A numeric vector of integer indices for the diagonal elements.
##' @keywords internal
vech_diag_indices <- function(d) {
	if (d == 0) return(integer(0))
	if (d == 1) return(1L)
  
	# The position of the k-th diagonal is 1 + the number of elements
	# in all preceding columns (d, d-1, ..., d-k+2).
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
   	if (d == 0) return(integer(0))
    	if (d == 1) return(1L) 

	list(
		i = unlist(lapply(1:d, function(j) j:d)),
		j = unlist(lapply(1:d, function(j) rep(j, d - j + 1)))
	)
}

##' @title Reconstruct Cholesky Factor from Parameters
##'
##' @description Internal helper function that takes an `UnstructuredCov` object
##' and reconstructs the dense lower-triangular Cholesky factor `L` from the
##' unconstrained parameter vector stored in the object's `@parameters` slot.
##'
##' @param object An S4 object inheriting from `VirtualCovariance`, which contains
##'  `@parameters` and `@dimension` slots.
##' @return A `dtrMatrix` object representing the Cholesky factor `L`.
##' @keywords internal
get_chol_from_params <- function(object) {
	d <- object@dimension
	if (d == 0) return(new("dtrMatrix", uplo = "L", Dim = c(0L, 0L)))

	L_vec <- object@parameters
	diag_indices <- vech_diag_indices(d)
	L_vec[diag_indices] <- exp(L_vec[diag_indices])

	L <- matrix(0, nrow = d, ncol = d)
	if (d == 1) {
		L[1,1] <- L_vec
	} else {
	vech_indices <- get_vech_indices(d)
	L[cbind(vech_indices$i, vech_indices$j)] <- L_vec
	}

	as(L, "dtrMatrix")
}

##' @title Get Interpretable Parameters for Compound Symmetry
##'
##' @description Internal helper function to back-transform the unconstrained
##' parameters of a `CompoundSymmetryCov` object into their interpretable
##' form (variance and correlation).
##'
##' @param object A `CompoundSymmetryCov` object.
##' @return A named list with two components: `variance` and `correlation`.
##' @keywords internal
get_cs_interpretable_params <- function(object) {
    list(
        variance = exp(object@parameters[1]),
        correlation = tanh(object@parameters[2])
    )
}
