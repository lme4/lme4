## Generate sparse multidimensional Gaussian quadrature grids --->  ../man/GQdk.Rd
## Unused currently; rather GHrule() --> ./GHrule.R
GQdk <- function(d=1L, k=1L) {
    stopifnot(0L < (d <- as.integer(d)[1]),
              d <= 20L,
              0L < (k <- as.integer(k)[1]),
	      k <= length(GQNd <- GQN[[d]]))## -> GQN, stored in ./sysdata.rda
    tmat <- t(GQNd[[k]])
    rperms <- lapply(.Call(allPerm_int, seq_len(d) + 1L), function(v) c(1L, v))
    dd <- unname(as.matrix(do.call(expand.grid, c(rep.int(list(c(-1,1)), d), KEEP.OUT.ATTRS=FALSE))))
    unname(unique(t(do.call(cbind,
                            lapply(as.data.frame(t(cbind(1, dd))),
                                   "*", e2=do.call(cbind, lapply(rperms, function(ind) tmat[ind,])))))))
}
