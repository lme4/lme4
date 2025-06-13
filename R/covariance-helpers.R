##' @title Get Vech Diagonal Indices
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
  if (d == 0) {
    return(integer(0))
  }
  # The position of the k-th diagonal is 1 + the number of elements
  # in all preceding columns (d, d-1, ..., d-k+2).
  1L + c(0, cumsum(d:(d - d + 2)))
}

##' @title Get Vech Indices
##' @description Internal helper to get row and column indices for a 'vech'
##' (column-major, lower-triangular) ordered vector.
##'
##' @param d An integer, the dimension of the square matrix.
##' @return A list with components 'i' (row indices) and 'j' (column indices).
##' @keywords internal
get_vech_indices <- function(d) {
  if (d == 0) {
    return(list(i = integer(0), j = integer(0)))
  }
  list(
    i = unlist(lapply(1:d, function(j) j:d)),
    j = unlist(lapply(1:d, function(j) rep(j, d - j + 1)))
  )
}
