## --> ../man/NelderMead.Rd
Nelder_Mead <- function(fn, par, lower=rep.int(-Inf, n),
                        upper=rep.int(Inf, n), control=list()) {
    n <- length(par)
    if (is.null(xst <- control[["xst"]])) xst <- rep.int(0.02,n)
    if (is.null(xt <- control[["xt"]])) xt <- xst*5e-4

    control[["xst"]] <- control[["xt"]] <- NULL

    ## mapping between simpler 'verbose' setting (0=no printing, 1=20, 2=10, 3=1)
    ##  and internal 'iprint' control (frequency of printing)
    if (is.null(verbose <- control[["verbose"]])) verbose <- 0
    control[["verbose"]] <- NULL
    if (is.null(control[["iprint"]])) {
      control[["iprint"]] <- switch(as.character(min(as.numeric(verbose),3L)),
                                    "0"=0, "1"=20,"2"=10,"3"=1)
    }
    stopifnot(is.function(fn),
              length(formals(fn)) == 1L,
              (n <- length(par <- as.numeric(par))) == length(lower <- as.numeric(lower)),
              length(upper <- as.numeric(upper)) == n,
              length(xst <- as.numeric(xst)) == n,
              all(xst != 0),
              length(xt <- as.numeric(xt)) == n)
    ## "NelderMead" reference class and constructor: --> ./AllClass.R :
    nM <- NelderMead$new(lower=lower, upper=upper, x0=par, xst=xst, xt=xt)
    cc <- do.call(function(iprint = 0L, maxfun = 10000L, FtolAbs = 1e-5,
                           FtolRel = 1e-15, XtolRel = 1e-7,
                           MinfMax= -.Machine$double.xmax, warnOnly=FALSE, ...) {
        if (length(list(...))>0) warning("unused control arguments ignored")
        list(iprint=iprint, maxfun=maxfun, FtolAbs=FtolAbs, FtolRel=FtolRel,
             XtolRel=XtolRel, MinfMax=MinfMax, warnOnly=warnOnly)
    }, control)
    nM$setFtolAbs(cc$FtolAbs)
    nM$setFtolRel(cc$FtolRel)
    nM$setIprint (cc$iprint)
    nM$setMaxeval(cc$maxfun)
    nM$setMinfMax(cc$MinfMax)
    it <- 0
    repeat {
        it <- it + 1
        nMres <- nM$newf(fn(nM$xeval()))
        if (nMres != 0L) break
    }

    cmsg <- "reached max evaluations"
    if (nMres == -4) {
        ## map max evals from error to warning
        cmsg <- warning(sprintf("failure to converge in %d evaluations",cc$maxfun))
        nMres <- 4
    }
                                                        ## nMres:
    msgvec <- c("nm_forced",                            ## -3
                "cannot generate a feasible simplex",   ## -2
                "initial x is not feasible",            ## -1
                "active",                               ## 0 (active)
                "objective function went below allowed minimum",           ## 1 (minf_max)
                "objective function values converged to within tolerance", ## 2 (fcvg)
                "parameter values converged to within tolerance",          ## 3 (xcvg)
                cmsg)

    if (nMres < 0) { ## i.e., in {-3, -2, -1}
        (if(cc$warnOnly) warning else stop)( msgvec[nMres+4] )
    }

    list(fval = nM$value(), par = nM$xpos(),
	 convergence = pmin(0, nMres), # positive nMres is also 'convergence'
	 NM.result = nMres, `message` = msgvec[nMres+4],
	 control = c(cc, xst=xst, xt=xt), feval = it)
}
