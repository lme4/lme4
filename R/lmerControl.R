namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}
## TESTING:
## a <- b <- c <- 1
## namedList(a,b,c)
## namedList(a,b,d=c)
## namedList(e=a,f=b,d=c)

##' By extracting "checking options" from \code{nms}, this function
##' implicitly defines what "checking options" are.
##'
##' @title Extract "checking options" or "checking arguments" (-> ../man/lmerControl.Rd)
##' @param nms character vector
##' @return those elements of \code{nms} which are "checking options"
##' @author Martin Maechler
.get.checkingOpts <- function(nms)
    nms[grepl("^check\\.(?!conv)", nms, perl=TRUE) | nms == "ensureXrank"]


## DOC: ../man/lmerControl.Rd
lmerControl <- function(optimizer="bobyqa",#instead of "Nelder_Mead",
                        restart_edge=TRUE,
                        ## don't call this "check." -- will fail
                        ## automatic check-option-checking in
                        ## inst/tests/test-lmer.R
                        boundary.tol=1e-5,
                        calc.derivs=TRUE,
                        use.last.params=FALSE,
                        sparseX=FALSE,
                        ## input checking options:
                        check.nobs.vs.rankZ="warningSmall",
                        check.nobs.vs.nlev="stop",
                        check.nlev.gtreq.5="ignore",
                        check.nlev.gtr.1="stop",
                        check.nobs.vs.nRE="stop",
                        ensureXrank=TRUE,
                        check.formula.LHS="stop",
                        check.conv.grad="warning",
                        check.conv.singular="ignore",
                        check.conv.hess="warning",
                        optCtrl = list()
                        )
{
    ## FIXME: is there a better idiom?  match.call() ?
    ## fill in values from options, but **only if not specified explicitly in arguments**
    ##  (ugh ... is there a better way to do this?  mapply() is clunky:
    ##  http://stackoverflow.com/questions/16276667/using-apply-with-assign-in-r
    stopifnot(is.list(optCtrl))

    if (!is.null(lmerOpts <- getOption("lmerControl"))) {
        for (arg in .get.checkingOpts(names(lmerOpts))) {
            if (do.call(missing,list(arg))) ## only if missing from explicit arguments
                assign(arg,lmerOpts[[arg]])
        }
    }
    structure(namedList(optimizer,
                        restart_edge,
                        boundary.tol,
                        calc.derivs,
                        use.last.params,
                        checkControl =
                        namedList(check.nobs.vs.rankZ,
                                  check.nobs.vs.nlev,
                                  check.nlev.gtreq.5,
                                  check.nlev.gtr.1,
                                  check.nobs.vs.nRE,
                                  ensureXrank,
                                  check.formula.LHS),
                        checkConv=
                        namedList(check.conv.grad,
                                  check.conv.singular,
                                  check.conv.hess),
                        optCtrl=optCtrl),
              class = c("lmerControl", "merControl"))
}

glmerControl <- function(optimizer=c("bobyqa","Nelder_Mead"),
                         restart_edge=FALSE,
                         boundary.tol=1e-5,
                         calc.derivs=TRUE,
                         use.last.params=FALSE,
                         sparseX=FALSE,
                         tolPwrss = 1e-7,
                         compDev = TRUE,
                         ## input checking options
                         check.nobs.vs.rankZ="warningSmall",
                         check.nobs.vs.nlev="stop",
                         check.nlev.gtreq.5="ignore",
                         check.nlev.gtr.1="stop",
                         check.nobs.vs.nRE="stop",
                         ## convergence checking options
                         check.conv.grad="warning",
                         check.conv.singular="ignore",
                         check.conv.hess="warning",
                         ensureXrank=TRUE,
                         check.formula.LHS="stop",
                         optCtrl = list())
{
    ## FIXME: should try to modularize/refactor/combine with lmerControl if possible
    ## but note different defaults
    ##                lmer        glmer
    ## optimizer    Nelder_Mead  c(Nelder_Mead,bobyqa)
    ## tolPwrss     N/A          1e-7
    ## compDev      N/A          TRUE
    ##
    ## (and possible future divergence)
    stopifnot(is.list(optCtrl))
    if (length(optimizer)==1) {
	optimizer <- replicate(2,optimizer) # works evevn when optimizer is function
    }
    if (!is.null(glmerOpts <- getOption("glmerControl"))) {
        for (arg in .get.checkingOpts(names(glmerOpts))) {
            if (do.call(missing,list(arg))) ## only if missing from explicit arguments
                assign(arg, glmerOpts[[arg]])
        }
    }
    if (use.last.params && calc.derivs)
        warning("using ",shQuote("use.last.params"),"=TRUE and ",
                shQuote("calc.derivs"),"=TRUE with ",shQuote("glmer"),
                " will not give backward-compatible results")
    structure(namedList(optimizer,
                        calc.derivs,
                        use.last.params,
			restart_edge,
                        boundary.tol,
			tolPwrss,
			compDev,
			checkControl=
			namedList(check.nobs.vs.rankZ,
                                  check.nobs.vs.nlev,
				  check.nlev.gtreq.5,
				  check.nlev.gtr.1,
                                  check.nobs.vs.nRE,
                                  ensureXrank,
                                  check.formula.LHS),
                        checkConv=
                        namedList(check.conv.grad,
                                  check.conv.singular,
                                  check.conv.hess),
                        optCtrl=optCtrl),
	      class = c("glmerControl", "merControl"))
}


##' @rdname lmerControl
##' @export
nlmerControl <- function(optimizer="Nelder_Mead",
                         tolPwrss = 1e-10,
                         optCtrl = list())
{
    stopifnot(is.list(optCtrl))
    if (length(optimizer)==1) {
        optimizer <- replicate(2,optimizer)
    }
    structure(namedList(optimizer,
			tolPwrss,
			optCtrl=optCtrl),
	      class = c("nlmerControl", "merControl"))
}
