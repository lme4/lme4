namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}

##' By extracting "checking options" from \code{nms}, this function
##' implicitly defines what "checking options" are.
##'
##' @title Extract "checking options" or "checking arguments" (-> ../man/lmerControl.Rd)
##' those starting with "check." but *not* the "check.conv.." nor "check.rankX.." ones:
##' @param nms character vector
##' @return those elements of \code{nms} which are "checking options"
##' @author Martin Maechler
.get.checkingOpts <- function(nms)
    nms[grepl("^check\\.(?!conv|rankX|scaleX)", nms, perl=TRUE)]


##' Check check.conv.*() options and produce good error message
chk.convOpt <- function(opt) {
    cnm <- deparse(nm <- substitute(opt))[[1]]
    if(!is.list(opt)) stop("check.conv* option ", cnm, " must be a list")
    if(!is.character(opt$action))
	stop("check.conv* option ", cnm, " has no 'action' string")
    if(!is.numeric(tol <- opt$tol))
	stop("check.conv* option ", cnm, " must have a numeric 'tol' component")
    if(length(tol) != 1 || tol < 0)
	stop("check.conv* option ", cnm, "$tol must be number >= 0")
    if(!is.null(relTol <- opt$relTol))
	stopifnot(is.numeric(relTol), length(relTol) == 1, relTol >= 0)
    invisible()
}

##' Exported constructor for the user calling  *lmerControl():
.makeCC <- function(action, tol, relTol, ...) {
    stopifnot(is.character(action), length(action) == 1)
    if(!is.numeric(tol))
	stop("must have a numeric 'tol' component")
    if(length(tol) != 1 || tol < 0)
	stop("'tol' must be number >= 0")
    mis.rt <- missing(relTol)
    if(!mis.rt && !is.null(relTol))
	stopifnot(is.numeric(relTol), length(relTol) == 1, relTol >= 0)
    ## and return the list,  the  "..." just being appended unchecked
    c(list(action = action, tol = tol),
      if(!mis.rt) list(relTol = relTol),
      list(...))
}

##' Internal utility :	Allow check.conv.* to be a string
chk.cconv <- function(copt, callingFn) {
    cnm <- deparse(substitute(copt))
    if(is.character(copt)) {
	def <- eval(formals(callingFn)[[cnm]])
	def$action <- copt
	assign(cnm, def, envir=sys.frame(sys.parent()))
    } else chk.convOpt(copt)
}


## DOC: ../man/lmerControl.Rd
lmerControl <-
    function(optimizer="bobyqa",#instead of "Nelder_Mead",
	     restart_edge=TRUE,
	     ## don't call this "check." -- will fail
	     ## automatic check-option-checking in
	     ## inst/tests/test-lmer.R
	     boundary.tol=1e-5,
	     calc.derivs=TRUE,
	     use.last.params=FALSE,
	     sparseX=FALSE,
	     ## input checking options:
	     check.nobs.vs.rankZ="ignore", ## "warningSmall",
	     check.nobs.vs.nlev="stop",
	     check.nlev.gtreq.5="ignore",
	     check.nlev.gtr.1="stop",
	     check.nobs.vs.nRE="stop",
	     check.rankX = c("message+drop.cols",
                             "silent.drop.cols", "warn+drop.cols",
                 	     "stop.deficient", "ignore"),
	     check.scaleX = c("warning","stop","silent.rescale",
                              "message+rescale","warn+rescale","ignore"),
	     check.formula.LHS = "stop",
	     ## convergence options
	     check.conv.grad	 = .makeCC("warning", tol = 2e-3, relTol = NULL),
	     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
	     check.conv.hess	 = .makeCC(action = "warning", tol = 1e-6),
	     optCtrl = list()
	     )
{
    ## FIXME: is there a better idiom?  match.call() ?
    ## fill in values from options, but **only if not specified explicitly in arguments**
    ##  (ugh ... is there a better way to do this?  mapply() is clunky:
    ##  http://stackoverflow.com/questions/16276667/using-apply-with-assign-in-r
    stopifnot(is.list(optCtrl))

    if (!is.null(lmerOpts <- getOption("lmerControl"))) {
        nn <- names(lmerOpts)
        nn.ok <- .get.checkingOpts(names(lmerOpts))
        if (length(nn.ignored <- setdiff(nn,nn.ok))>0) {
            warning("some options in ",shQuote("getOption('lmerControl')"),
                    " ignored : ",paste(nn.ignored,collapse=", "))
        }
        for (arg in nn.ok) {
            if (do.call(missing,list(arg))) ## only if missing from explicit arguments
                assign(arg,lmerOpts[[arg]])
        }
    }
    check.rankX <- match.arg(check.rankX)# ==> can abbreviate
    check.scaleX <- match.arg(check.scaleX)# ==> can abbreviate

    ## compatibility and convenience, caller can specify action string only:
    me <- sys.function()
    chk.cconv(check.conv.grad,	   me)
    chk.cconv(check.conv.singular, me)
    chk.cconv(check.conv.hess	 , me)

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
				  check.rankX,
                                  check.scaleX,
                                  check.formula.LHS),
                        checkConv=
                        namedList(check.conv.grad,
                                  check.conv.singular,
                                  check.conv.hess),
                        optCtrl=optCtrl),
              class = c("lmerControl", "merControl"))
}

glmerControl <-
    function(optimizer=c("bobyqa","Nelder_Mead"),
	     restart_edge=FALSE,
	     boundary.tol=1e-5,
	     calc.derivs=TRUE,
	     use.last.params=FALSE,
	     sparseX=FALSE,
	     tolPwrss = 1e-7,
	     compDev = TRUE,
	     ## input checking options
	     check.nobs.vs.rankZ="ignore", ## "warningSmall",
	     check.nobs.vs.nlev="stop",
	     check.nlev.gtreq.5="ignore",
	     check.nlev.gtr.1="stop",
	     check.nobs.vs.nRE="stop",
	     check.rankX = c("message+drop.cols",
	     "silent.drop.cols", "warn+drop.cols",
	     "stop.deficient", "ignore"),
	     check.scaleX = "warning",
	     check.formula.LHS="stop",
	     ## convergence checking options
	     check.conv.grad	 = .makeCC("warning", tol = 1e-3, relTol = NULL),
	     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
	     check.conv.hess	 = .makeCC(action = "warning", tol = 1e-6),
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
        nn <- names(glmerOpts)
        nn.ok <- .get.checkingOpts(names(glmerOpts))
        if (length(nn.ignored <- setdiff(nn,nn.ok))>0) {
            warning("some options in ",shQuote("getOption('glmerControl')"),
                    " ignored : ",paste(nn.ignored,collapse=", "))
        }
        for (arg in nn.ok) {
            if (do.call(missing,list(arg))) ## only if missing from explicit arguments
                assign(arg, glmerOpts[[arg]])
        }
    }
    check.rankX <- match.arg(check.rankX)# ==> can abbreviate

    ## compatibility and convenience, caller can specify action string only:
    me <- sys.function()
    chk.cconv(check.conv.grad,	   me)
    chk.cconv(check.conv.singular, me)
    chk.cconv(check.conv.hess	 , me)

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
				  check.rankX,
                                  check.scaleX,
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
