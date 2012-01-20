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
