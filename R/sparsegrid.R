##' Reproduce the sparse grid for Gauss-Hermite quadrature from stored values
##'
##' The GQN list of lists stores the condensed arrays from \url{www.sparse-grids.de}.
##' This function recreates the full array by applying permutations and +/- multipliers.
##' @title Multidimensional Gauss-Hermite sparse interpolation grid
##' @param d dimension of grid
##' @param k order of approximation.  The grid will provide accurate integrals for
##'        polynomials of total order <= 2k - 1 multiplied by the kernel
##' @return the grid as a matrix.  The first column contains the weights.  The
##'        remaining d columns are the coordinates on the z scale.
##' @author Douglas Bates
##' @examples
##' GQdk(3, 5)
##' @export
GQdk <- function(d=1L, k=1L) {
    stopifnot(0L < (d <- as.integer(d)[1]),
              d <= 20L,
              0L < (k <- as.integer(k)[1]),
              k <= length(GQNd <- GQN[[d]]))
    tmat   <- t(GQNd[[k]])
    dseq   <- seq_len(d)
    rperms <- lapply(as.data.frame(.Call(allPerm_int, dseq + 1L)), function(v) c(1L, v))
    unname(unique(t(do.call(cbind,
                            lapply(as.data.frame(t(cbind(1,
                                                         as.matrix(do.call(expand.grid,
                                                                           lapply(dseq,
                                                                                  function(i) c(-1,1))))))),
                                   "*", e2=do.call(cbind, lapply(rperms, function(ind) tmat[ind,])))))))
}
