library("lme4")
library("numDeriv")
namedList <- lme4:::namedList ## from lmerControl.R

simfun <- function(size=nrow(data),seed=NULL,
                   data=Penicillin,
                   formula=diameter ~ 1 + (1|plate) + (1|sample),
                   fn = lmer,
                   ...) {
    if (!is.null(seed)) set.seed(seed)
    ## bootstrap or pseudo-bootstrap
    sdat <- data[sample(x=nrow(data),size=size,replace=TRUE),]
    t1 <- system.time(fm2 <- fn(formula, sdat, ...))
    derivs <- fm2@optinfo$derivs
    ## copied from lme4:::checkConv
    mineig <- min(eigen(fm2@optinfo$derivs$Hessian,only.values=TRUE)$value)
    dd <- as.function(fm2)  ## FIXME for glmer()
    t2 <- system.time(hh <- hessian(dd,getME(fm2,"theta")))
    mineigND <- min(eigen(hh,only.values=TRUE)$value)
    scgrad <- tryCatch(with(derivs, solve(Hessian, gradient)), 
                       error = function(e) e)
    if (inherits(scgrad, "error")) {
        wstr <- "unable to evaluate scaled gradient"
        maxscgrad <- maxmingrad <- NA
    } else {
        maxscgrad <- max(abs(scgrad))
        mingrad <- pmin(abs(scgrad), abs(derivs$gradient))
        maxmingrad <- max(mingrad)
    }
    th <- getME(fm2,"theta")
    names(th) <- paste0("theta",seq_along(th))
    res <- unlist(namedList(maxgrad=max(abs(derivs$gradient)),
                            maxscgrad,
                            maxmingrad,
                            mineig,
                            mineigND,
                            th,
                            t.tot=unname(t1["elapsed"]),
                            t.hessian=unname(t2["elapsed"])))
    return(res)
}
