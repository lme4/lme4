##' Generate the sparse multidimensional Gaussian quadrature grids
##'
##' @title Sparse Gaussian Quadrature grid
##' @param d integer scalar - the dimension of the function to be
##'    integrated with respect to the standard \code{d}-dimensional
##'    Gaussian density
##' @param k integer scalar - the order of the grid.  A grid of order
##'    \code{k} provides an exact result for a polynomial of total order
##'    of \code{2k - 1} or less multiplied by the 
##' @return a matrix with \code{d + 1} columns.  The first column is
##'    the weights and the remaining \code{d} columns are the node
##'    coordinates.
##' @note The number of nodes gets very large very quickly with
##'    increasing \code{d} and \code{k}.  See the charts at
##'    \url{http://www.sparse-grids.de}.
##' @examples
##' GQdk(2,5)
##' @export
GQdk <- function(d=1L, k=1L) {
    stopifnot(0L < (d <- as.integer(d)[1]),
              d <= 20L,
              0L < (k <- as.integer(k)[1]),
              k <= length(GQNd <- GQN[[d]]))
    tmat   <- t(GQNd[[k]])
    dseq   <- seq_len(d)
    rperms <- lapply(.Call(allPerm_int, dseq + 1L), function(v) c(1L, v))
    unname(unique(t(do.call(cbind,
                            lapply(as.data.frame(t(cbind(1,
                                                         as.matrix(do.call(expand.grid,
                                                                           lapply(dseq,
                                                                                  function(i) c(-1,1))))))),
                                   "*", e2=do.call(cbind, lapply(rperms, function(ind) tmat[ind,])))))))
}
