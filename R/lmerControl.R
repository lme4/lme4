namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
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


merControl <-
    function(optimizer="bobyqa",
	     restart_edge=TRUE,
	     ## don't call this "check." -- will fail
	     ## automatic check-option-checking in
	     ## inst/tests/test-lmer.R
	     boundary.tol=1e-5,
	     calc.derivs=TRUE,
	     use.last.params=FALSE,
	     sparseX=FALSE,
             tolPwrss = 1e-7,       ## GLMM only
             compDev = TRUE,        ## GLMM only
             nAGQ0initStep = TRUE,  ## GLMM only
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
             check.response.not.const = "stop",  ## GLMM only
	     ## convergence options
	     check.conv.grad	 = .makeCC("warning", tol = 0.01, relTol = NULL),
	     check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
	     check.conv.hess	 = .makeCC(action = "warning", tol = 1e-6),
	     optCtrl = list(),
             ## specify lmer vs glmer: LAST to avoid screwing up
             ##   positional matching
             type="lmer"  
             )

{

    ## FIXME: is there a better idiom?  match.call() ?
    ## fill in values from options, but **only if not specified explicitly in arguments**
    ##  (ugh ... is there a better way to do this?  mapply() is clunky:
    ##  http://stackoverflow.com/questions/16276667/using-apply-with-assign-in-r
    stopifnot(is.list(optCtrl))

    merControl <- paste0(type,"Control")
    if (!is.null(merOpts <- getOption(merControl))) {
        nn <- names(merOpts)
        nn.ok <- .get.checkingOpts(names(merOpts))
        if (length(nn.ignored <- setdiff(nn,nn.ok))>0) {
            warning("some options in ",shQuote(paste0("getOption('",merControl,"')")),
                    " ignored : ",paste(nn.ignored,collapse=", "))
        }
        for (arg in nn.ok) {
            if (do.call(missing,list(arg))) ## only if missing from explicit arguments
                assign(arg,merOpts[[arg]])
        }
    }
    check.rankX <- match.arg(check.rankX)# ==> can abbreviate
    check.scaleX <- match.arg(check.scaleX)# ==> can abbreviate

    ## compatibility and convenience, caller can specify action string only:
    me <- sys.function()
    chk.cconv(check.conv.grad,	   me)
    chk.cconv(check.conv.singular, me)
    chk.cconv(check.conv.hess	 , me)

    if (type=="glmer" && use.last.params && calc.derivs)
        warning("using ",shQuote("use.last.params"),"=TRUE and ",
                shQuote("calc.derivs"),"=TRUE with ",shQuote("glmer"),
                " will not give backward-compatible results")

    res <- namedList(optimizer,
                     calc.derivs,
                     use.last.params,
                     restart_edge,
                     boundary.tol)

    checkControl <- namedList(check.nobs.vs.rankZ,
                              check.nobs.vs.nlev,
                              check.nlev.gtreq.5,
                              check.nlev.gtr.1,
                              check.nobs.vs.nRE,
                              check.rankX,
                              check.scaleX,
                              check.formula.LHS)
    
    if (type=="glmer") {
        if (length(res$optimizer)==1) {
            res$optimizer <- replicate(2,res$optimizer) # works even when optimizer is function
        }

        res <- c(res, namedList(tolPwrss,
                                compDev,
                                nAGQ0initStep))

        checkControl <- c(checkControl, namedList(check.response.not.const))
    }

    res <- c(res, namedList(checkControl,
                            checkConv=
                                namedList(check.conv.grad,
                                          check.conv.singular,
                                          check.conv.hess),
                            optCtrl))

    structure(res,class = c(merControl, "merControl"))
}


## DOC: ../man/lmerControl.Rd
lmerControl <- function(...) {
    mc <- match.call()
    ## FIXME: we need eval.parent(), but merControl isn't currently
    ## exported.  If we export it we need to document it (ugh).
    ## ::: triggers NOTEs
    ##  defined
    mc[[1]] <- quote(merControl)
    ## eval.parent(mc)
    eval.parent(mc)
}

## hack formals so that lmerControl matches documentation
## (at present we don't export merControl())
ff <- formals(merControl)
## drop glmer-specific controls and "type"
ff <- ff[!names(ff) %in% c("tolPwrss","compDev","nAGQ0initStep","check.response.not.const","type")]
formals(lmerControl) <- ff

## almost the same body but need to set 'type="glmer"' in arguments
glmerControl <- function(...) {
    mc <- match.call()
    mc[[1]] <- quote(merControl)
    mc[["type"]] <- "glmer"
    ## eval.parent(mc)
    eval.parent(mc)
}

ff <- formals(merControl)
## drop "type" (no LMM-specific controls at present)
ff <- ff[!names(ff) %in% "type"]

formals(glmerControl) <- ff

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
