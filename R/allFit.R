meth.tab.0 <- cbind(optimizer=
                      rep(c("bobyqa",
                            "Nelder_Mead",
                            "nlminbwrap",
                            "nmkbw",
                            "optimx",
                            "nloptwrap" ),
                          c(rep(1,5),2)),
                    method= c(rep("",4), "L-BFGS-B",
                            "NLOPT_LN_NELDERMEAD", "NLOPT_LN_BOBYQA"))

## ugh: hardcoded list (incomplete?) of allowable *control* options by optimizer
## could make more of an effort to match max-iterations/evals,
## (x|f)*(abs|rel) tolerance, ...
opt.ctrls <- list(bobyqa=c("npt","rhobeg","rhoend","iprint","maxfun"),
                  Nelder_Mead=c("iprint","maxfun","FtolAbs",
                                "FtolRel","XtolRel","MinfMax",
                                "xst","xt","verbose","warnOnly"),
                  nlminbwrap=c("eval.max","iter.max","trace","abs.tol",
                               "rel.tol","x.tol","xf.tol","step.min",
                               "step.max","sing.tol","scale.init",
                               "diff.g"),
                  optimx=c("trace","fnscale","parscale","ndeps",
                          "maxit","abstol","reltol","method"),
                  nloptwrap=c("minf_max","ftol_rel","ftol_abs",
                              "xtol_rel", "xtol_abs", "maxeval", "maxtime",
                              "algorithm"),
                  nmkbw=c("tol","maxfeval","regsimp","maximize",
                          "restarts.max","trace","maxfun"))

nmkbw <- function(fn,par,lower,upper,control) {
    if (length(par)==1) {
        res <- optim(fn=fn,par=par,lower=lower,upper=100*par,
                     method="Brent")
    } else {
        if (!is.null(control$maxfun)) {
            control$maxfeval <- control$maxfun
            control$maxfun <- NULL
        }
        res <- dfoptim::nmkb(fn=fn,par=par,
                             lower=lower,upper=upper,control=control)
    }
    res$fval <- res$value
    res
}


##' Attempt to re-fit a [g]lmer model with a range of optimizers.
##' The default is to use all known optimizers for R that satisfy the
##' requirements (do not require explicit gradients, allow
##' box constraints), in three categories; (i) built-in
##' (minqa::bobyqa, lme4::Nelder_Mead, nlminbwrap), (ii) wrapped via optimx
##' (most of optimx's optimizers that allow box constraints require
##' an explicit gradient function to be specified; the two provided
##' here are really base R functions that can be accessed via optimx,
##' (iii) wrapped via nloptr; (iv)
##'
##' @param object a fitted model
##' @param meth.tab a matrix (or data.frame) with columns
##' - method  the name of a specific optimization method to pass to the optimizer
##'           (leave blank for built-in optimizers)
##' - optimizer  the \code{optimizer} function to use
##' @param verbose print progress messages?
##' @return a list of fitted \code{merMod} objects
##' @seealso slice, slice2D in the bbmle package
##' @examples
##' library(lme4)
##' gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'                  data = cbpp, family = binomial)
##' gm_all <- allFit(gm1, parallel=TRUE)
##' ss <- summary(gm_all)
##' ss$fixef               ## extract fixed effects
##' ss$llik                ## log-likelihoods
##' ss$sdcor               ## SDs and correlations
##' ss$theta               ## Cholesky factors
##' ss$which.OK            ## which fits worked

allFit <- function(object, meth.tab = NULL,
                   data=NULL,
                   verbose=TRUE,
                   show.meth.tab = FALSE,
                   maxfun=1e5,
                   parallel = c("no", "multicore", "snow"),
                   ncpus = getOption("allFit.ncpus", 1L),
                   cl = NULL) {

    if (is.null(meth.tab)) {
        meth.tab <- meth.tab.0
    }
    if (!requireNamespace("dfoptim")) {
        meth.tab <- meth.tab[meth.tab.0[,"optimizer"] != "nmkbw",]
    }
    if (!requireNamespace("optimx")) {
        meth.tab <- meth.tab[meth.tab.0[,"optimizer"] != "optimx",]
    }
    if (show.meth.tab) {
        return(meth.tab)
    }
    stopifnot(length(dm <- dim(meth.tab)) == 2, dm[1] >= 1, dm[2] >= 2,
	      is.character(optimizer <- meth.tab[,"optimizer"]),
	      is.character(method    <- meth.tab[,"method"]))

    parallel <- match.arg(parallel)
    do_parallel <- have_mc <- have_snow <- NULL # "-Wall"
    eval(initialize.parallel) # (parallel, ncpus)   --> ./utilities.R
    ## |--> (do_parallel, have_mc, have_snow)

    fit.names <- gsub("\\.$", "", paste(optimizer, method, sep="."))
    ffun <- local({
        ## required local vars
        object
        verbose
        fit.names
        optimizer
        method
        maxfun
        function(..i) {
            if (verbose) cat(fit.names[..i],": ")
            ctrl <- getCall(object)$control
            ## NB:  'ctrl' must become a correct *argument* list for g?lmerControl()
            if(is.null(ctrl)) {
                ctrl <- list(optimizer=optimizer[..i])
            } else {
                if(is.call(ctrl)) # typically true
                    ctrl <- lapply(as.list(ctrl)[-1], eval)
                ctrl$optimizer <- optimizer[..i]
            }
            ## add method/algorithm to optCtrl as necessary
            mkOptCtrl <- function(...) {
                x <- list(...)
                cc <- ctrl$optCtrl
                if (is.null(cc)) return(x)  ## easy! no controls specified
                for (n in names(x)) { ## replace existing values
                    cc[[n]] <- x[[n]]
                }
                cc
            }
            sanitize <- function(x,okvals) {
                if (is.null(x)) return(NULL)
                if (is.null(okvals)) return(x)
                x <- x[intersect(names(x),okvals)]
                ## ?? having names(control) be character(0)
                ##  screws up nmkbw ... ??
                if (length(names(x))==0) names(x) <- NULL
                x
            }
            ctrl$optCtrl <- switch(optimizer[..i],
                                   optimx    = mkOptCtrl(method   = method[..i]),
                                   nloptwrap = mkOptCtrl(algorithm= method[..i]),
                                   mkOptCtrl("maxfun"=maxfun))
            ctrl$optCtrl <- sanitize(ctrl$optCtrl,
                                     opt.ctrls[[optimizer[..i]]])
            ctrl <- do.call(if(isGLMM(object)) glmerControl else lmerControl, ctrl)
            tt <- system.time(rr <- tryCatch(update(object, control = ctrl),
                                             error = function(e) e))
            attr(rr, "optCtrl") <- ctrl$optCtrl # contains crucial info here
            attr(rr, "time") <- tt  # store timing info
            if (verbose) cat("[OK]\n")
            return(rr)
        }
    })

    seq_fit <- seq_along(fit.names)
    res <- if (do_parallel) {
               if (have_mc) {
                   parallel::mclapply(seq_fit,
                                      ffun, mc.cores = ncpus)
               } else if(have_snow) {
                   if(is.null(cl)) {
                       cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                       ## consider exporting data/package ?
                       ## parallel::clusterEvalQ(cl,library("lme4"))
                       ## try to get data and export it?
                       ## parallel::clusterExport(cl,??)
                       res <- parallel::parLapply(cl, seq_fit, ffun)
                       parallel::stopCluster(cl)
                       res
                   } else parallel::parLapply(cl, seq_fit, ffun)
               } else {
                   warning("'do_parallel' is true, but 'have_mc' and 'have_snow' are not.  Should not happen!")
                   ## or stop()  or  we could silently use lapply(..)
                   setNames(as.list(fit.names), fit.names)
               }
           } else
               lapply(seq_fit, ffun)

    names(res) <- fit.names
    structure(res, class = "allFit", fit = object, sessionInfo =  sessionInfo(),
              data = data # is dropped if NULL
              )
}

print.allFit <- function(x, width=80, ...) {
    cat("original model:\n")
    f <- attr(x,"fit")
    ss <- function(x) {
        if (nchar(x)>width) {
            strwrap(paste0(substr(x,1,width-3),"..."))
        } else x
    }
    ff <- ss(safeDeparse(formula(f)))
    cat(ff,"\n")
    cat("data: ",deparse(getCall(f)$data),"\n")
    cat("optimizers (",length(x),"): ",
        ss(paste(names(x),collapse=", ")),"\n",
        sep="")
    which.bad <- vapply(x,FUN=is,"error",FUN.VALUE=logical(1))
    if ((nbad <- sum(which.bad))>0) {
        cat(nbad,"optimizer(s) failed\n")
    }
    cat("differences in negative log-likelihoods:\n")
    nllvec <- -vapply(x[!which.bad],logLik,numeric(1))
    cat("max=",signif(max(nllvec-min(nllvec)),3),
        "; std dev=",signif(sd(nllvec),3), "\n")
    ## FIXME: print magnitudes of parameter diffs
    ## cat("differences in parameters:\n")
    ## ss <- summary(x)
    ## allpars <- cbind(ss$fixef, ss$sdcor)
    ## par_max <-
    invisible(x)
}

summary.allFit <- function(object, ...) {
    afun <- function(x, FUN, ...) {
        f1 <- FUN(x[[1]], ...)
        nm <- names(f1)
        n <- length(f1)
        res <- vapply(x, FUN, numeric(n), ...)
        if (!is.null(dim(res))) {
            res <- t(res)
        } else {
            res <- as.matrix(res)
            colnames(res) <- nm
        }
        res
    }
    which.OK <- !vapply(object, is, "error", FUN.VALUE=logical(1))
    objOK <- object[which.OK]
    msgs <- lapply(objOK, function(x) x@optinfo$conv$lme4$messages)
    fixef <- afun(objOK, fixef)
    llik <- vapply(objOK, logLik, numeric(1))
    times <- afun(objOK, attr, "time")
    feval <- vapply(objOK, function(x) x@optinfo$feval, numeric(1))
    vfun <- function(x) as.data.frame(VarCorr(x))[["sdcor"]]
    sdcor <- afun(objOK, vfun)
    theta <- afun(objOK, getME, name="theta")
    cnm <- tnames(objOK[[1]])
    if (sigma(objOK[[1]])!=1) cnm <- c(cnm,"sigma")
    colnames(sdcor) <- unname(cnm)
    sdcor <- as.data.frame(sdcor)
    res <- namedList(which.OK, msgs, fixef, llik, sdcor, theta, times, feval)
    class(res) <- "summary.allFit"
    res
}

## should add a print method for summary:
##  * fixed effects, random effects: summary of differences?

## not yet ...
plot.allFit <- function(x, abbr=16, ...) {
    values <- opt <- NULL ## R code check/non-standard evaluation
     if (! (requireNamespace("ggplot2"))) {
         stop("ggplot2 package must be installed to plot allFit objects")
     }
     aes <- NULL ## code check
     ss <- summary(x)
     ff <- stack(as.data.frame(ss$fixef))
     ff$opt <- rep(rownames(ss$fixef),length.out=nrow(ff))
     if (!is.null(abbr)) ff$opt <- abbreviate(ff$opt, minlength=abbr)
     (ggplot2::ggplot(ff, aes(values, opt, colour=opt))
         + ggplot2::geom_point()
         + ggplot2::facet_wrap(~ind,scale="free")
         + ggplot2::theme(legend.position="none")
     )
}
