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
##' @param m a fitted model
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

allFit <- function(m, meth.tab = NULL,
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
    if (show.meth.tab) {
        return(meth.tab)
    }

    do_parallel <- have_mc <- have_snow <- NULL
    eval(initialize.parallel)

    stopifnot(length(dm <- dim(meth.tab)) == 2, dm[1] >= 1, dm[2] >= 2,
	      is.character(optimizer <- meth.tab[,"optimizer"]),
	      is.character(method    <- meth.tab[,"method"]))
    fit.names <- gsub("\\.$","",paste(optimizer, method, sep="."))
    res <- setNames(as.list(fit.names), fit.names)
    ffun <- local({
        ## required local vars
        m
        verbose
        fit.names
        optimizer
        method
        maxfun
        function(i) {
            if (verbose) cat(fit.names[i],": ")
            ctrl <- getCall(m)$control
            if (is.null(ctrl)) {
                ctrl <- list(optimizer=optimizer[i])
            } else {
                ctrl$optimizer <- optimizer[i]
            }
            ctrl$optCtrl <- switch(optimizer[i],
                                   optimx    = list(method   = method[i]),
                                   nloptWrap = list(algorithm= method[i]),
                                   list(maxfun=maxfun))
            ctrl <- do.call(if(isGLMM(m)) glmerControl else lmerControl, ctrl)
            tt <- system.time(rr <- tryCatch(update(m, control = ctrl),
                                             error = function(e) e))
            attr(rr, "optCtrl") <- ctrl$optCtrl # contains crucial info here
            attr(rr, "time") <- tt  # store timing info
            if (verbose) cat("[OK]\n")
            return(rr)
        }
    })

    res <- if (do_parallel) {
               if (have_mc) {
                   parallel::mclapply(seq_along(fit.names),
                                                ffun, mc.cores = ncpus)
        } else if (have_snow) {
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                res <- parallel::parLapply(cl, seq_along(fit.names), ffun)
                parallel::stopCluster(cl)
                res
            } else parallel::parLapply(cl, seq_along(fit.names), ffun)
        }
    } else lapply(seq_along(fit.names), ffun)

    structure(res, class = "allfit", fit = m, sessionInfo =  sessionInfo(),
              data = data # is dropped if NULL
              )
}

print.allfit <- function(x, width=80, ...) {
    cat("original model:\n")
    f <- attr(x,"fit")
    ss <- function(x) {
        if (nchar(x)>width) {
            strwrap(paste0(substr(x,1,width-3),"..."))
        } else x
    }
    ff <- ss(deparse(formula(f)))
    cat(ff,"\n")
    cat("data: ",deparse(getCall(f)$data),"\n")
    cat("optimizers (",length(x),"): ",
        ss(paste(names(x),collapse=", ")),"\n",
        sep="")
    which.bad <- sapply(x,is,"error")
    if ((nbad <- sum(which.bad))>0) {
        cat(nbad,"optimizers failed\n")
    }
    cat("differences in negative log-likelihoods:\n")
    nllvec <- -sapply(x[!which.bad],logLik)
    print(summary(nllvec-min(nllvec)))
}

summary.allfit <- function(object, ...) {
    which.OK <- !sapply(object, is, "error")
    objOK <- object[which.OK]
    msgs <- lapply(objOK, function(x) x@optinfo$conv$lme4$messages)
    fixef <- do.call(rbind, lapply(objOK, fixef))
    llik <- sapply(objOK, logLik)
    times <- t(sapply(objOK, attr, "time"))
    feval <- sapply(objOK, function(x) x@optinfo$feval)
    sdcor <- do.call(rbind, lapply(objOK, function(x)
        as.data.frame(VarCorr(x))[, "sdcor"]))
    theta <- do.call(rbind, lapply(objOK,
                                   function(x) getME(x, "theta")))
    cnm <- tnames(objOK[[1]])
    if (sigma(object[[1]])!=1) cnm <- c(cnm,"sigma")
    colnames(sdcor) <- unname(cnm)
    sdcor <- as.data.frame(sdcor)
    res <- namedList(which.OK, msgs, fixef, llik, sdcor, theta, times, feval)
    class(res) <- "summary.allFit"
    res
}

## should add a print method for summary:
##  * fixed effects, random effects: summary of differences?
