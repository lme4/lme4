if(FALSE) ## C++ code in ../src/mcmcsamp.cpp -- is also  #ifdef 0
# @S3method mcmcsamp merMod
mcmcsamp.merMod <- function(object, n=1L, verbose=FALSE, saveb=FALSE, ...) {
    n <- max(1L, as.integer(n)[1])
    dd <- getME(object, "devcomp")$dims
    ranef <- matrix(numeric(0), nrow = dd[["q"]], ncol = 0)
    if (saveb) ranef <- matrix(, nrow = dd[["q"]], ncol = n)
    sigma <- matrix(unname(sigma(object)), nrow = 1,
                    ncol = (if (dd[["useSc"]]) n else 0))
    ff <- fixef(object)
    fixef <- matrix(ff, nrow=dd[["p"]], ncol=n)
    rownames(fixef) <- names(ff)
    ## FIXME create a copy of the resp and pred modules
    ans <- new("merMCMC",
               Gp = object@Gp,
 #              ST = matrix(.Call(mer_ST_getPars, object), dd[["np"]], n),
               call = object@call,
               dims = object@dims,
               deviance = rep.int(unname(object@deviance[["ML"]]), n),
               fixef = fixef,
               nc = sapply(object@ST, nrow),
               ranef = ranef,
               sigma = sigma)
    .Call(mer_MCMCsamp, ans, object)
}
