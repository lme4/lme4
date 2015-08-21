## --> ../man/profile-methods.Rd

profnames <- function(object,signames) {
    useSc <- as.logical(object@devcomp$dims[["useSc"]])
    ntp <- length(getME(object,"theta"))
    nn <- if (signames) {
        sprintf(".sig%02d",seq(ntp))
    } else {
        tnames(object,old=FALSE,prefix=c("sd","cor"))
    }
    if (useSc) {
        nn <- c(nn,if (signames) ".sigma" else "sigma")
    }
    return(nn)
}

##' @importFrom splines backSpline interpSpline periodicSpline
##' @importFrom stats profile
##' @method profile merMod
##' @export
profile.merMod <- function(fitted,
                           which=NULL,
                           alphamax = 0.01,
			   maxpts = 100,
                           delta = NULL,
                           delta.cutoff = 1/8,
                           verbose=0, devtol=1e-9,
                           maxmult = 10,
                           startmethod = "prev",
                           optimizer = NULL,
                           control = NULL,
                           signames = TRUE,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("profile.ncpus", 1L), cl = NULL, ...)
{

    ## FIXME: allow choice of nextstep/nextstart algorithm?
    ## FIXME: allow selection of individual variables to profile by name?
    ## FIXME: allow for failure of bounds (non-pos-definite correlation matrices) when >1 cor parameter

    if (missing(parallel)) parallel <- getOption("profile.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    do_parallel <- (parallel != "no" && ncpus > 1L)
    if (do_parallel) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
	if (!(have_mc || have_snow))
	    do_parallel <- FALSE # (only for "windows")
    }

    if (is.null(optimizer)) optimizer <- fitted@optinfo$optimizer
    ## hack: doesn't work to set bobyqa parameters to *ending* values stored
    ## in @optinfo$control
    ignore.pars <- c("xst", "xt")
    control.internal <- fitted@optinfo$control
    if (length(ign <- which(names(control.internal) %in% ignore.pars)) > 0)
        control.internal <- control.internal[-ign]
    if (!is.null(control)) {
        i <- names(control)
        control.internal[[i]] <- control[[i]]
    }
    control <- control.internal
    ## parallel stuff copied from bootMer ...
    if (missing(parallel)) parallel <- getOption("profile.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncpus <- 1L
    }
    useSc <- isLMM(fitted) || isNLMM(fitted)
    dd <- devfun2(fitted,useSc,signames)
    ## FIXME: figure out to what do here ...
    if (isGLMM(fitted) && fitted@devcomp$dims[["useSc"]])
        stop("can't (yet) profile GLMMs with non-fixed scale parameters")
    stopifnot(devtol >= 0)
    base <- attr(dd, "basedev")
    thopt <- attr(dd, "thopt")
    stderr <- attr(dd, "stderr")
    pp <- environment(dd)$pp
    X.orig <- pp$X
    n <- environment(dd)$n
    p <- length(pp$beta0)

    opt <- attr(dd, "optimum")
    nptot <- length(opt)
    nvp <- nptot - p    # number of variance-covariance pars
    wi.vp <- seq_len(nvp)
    if(nvp > 0) fe.orig <- opt[- wi.vp]

    which <- get.which(which, nvp, nptot, names(opt), verbose)

    res <- c(.zeta = 0, opt)
    res <- matrix(res, nrow = maxpts, ncol = length(res),
                  dimnames = list(NULL, names(res)), byrow = TRUE)

    ## FIXME: why is cutoff based on nptot
    ## (i.e. boundary of simultaneous LRT conf region for nptot values)
    ##  when we are computing (at most) 2-parameter profiles here?
    cutoff <- sqrt(qchisq(1 - alphamax, nptot))

    if (is.null(delta))
        delta <- cutoff*delta.cutoff

    ## helper functions

    ## nextpar calculates the next value of the parameter being
    ## profiled based on the desired step in the profile zeta
    ## (absstep) and the values of zeta and column cc for rows
    ## r-1 and r.  The parameter may not be below lower (or above upper)
    nextpar <- function(mat, cc, r, absstep,
                        lower = -Inf, upper = Inf, minstep=1e-6) {
        rows <- r - (1:0)         # previous two row numbers
        pvals <- mat[rows, cc]
        zeta <- mat[rows, ".zeta"]
        num <- diff(pvals)
        if (is.na(denom <- diff(zeta)) || denom==0) {
            warning("Last two rows have identical or NA .zeta values: using minstep")
            step <- minstep
        } else {
            step <- absstep*num/denom
            if (step<0) {
                warning("unexpected decrease in profile: using minstep")
                step <- minstep
            } else {
                if (r>1) {
                    if (abs(step) > (maxstep <- abs(maxmult*num))) {
                        maxstep <- sign(step)*maxstep
                        if (verbose) cat(sprintf("capped step at %1.2f (multiplier=%1.2f > %1.2f)\n",
                                                 maxstep,abs(step/num),maxmult))
                        step <- maxstep
                    }
                }
            }
        }
        min(upper, max(lower, pvals[2] + sign(num) * step))
    }

    nextstart <- function(mat, pind, r, method) {
        ## FIXME: indexing may need to be checked (e.g. for fixed-effect parameters)
        switch(method,
               global= opt[seqpar1][-pind],  ## address opt, no zeta column
               prev  = mat[r,1+seqpar1][-pind],
               extrap = stop("not yet implemented"),## do something with mat[r-(1:0),1+seqnvp])[-pind]
               stop("invalid nextstart method"))
    }

    ## mkpar generates the parameter vector of theta and
    ## sigma from the values being profiled in position w
    mkpar <- function(np, w, pw, pmw) {
        par <- numeric(np)
        par[w] <- pw
        par[-w] <- pmw
        par
    }

    ## fillmat fills the third and subsequent rows of the matrix
    ## using nextpar and zeta
    ## FIXME:  add code to evaluate more rows near the minimum if that
    ##        constraint was active.
    fillmat <- function(mat, lowcut, upcut, zetafun, cc) {
        nr <- nrow(mat)
        i <- 2L
        while (i < nr && mat[i, cc] > lowcut && mat[i,cc] < upcut &&
                   (is.na(curzeta <- abs(mat[i, ".zeta"])) || curzeta <= cutoff)) {
            np <- nextpar(mat, cc, i, delta, lowcut, upcut)
            ns <- nextstart(mat, pind = cc-1, r=i, method=startmethod)
            mat[i + 1L, ] <- zetafun(np,ns)
            if (verbose>0) {
                cat(i,cc,mat[i+1L,],"\n")
            }
            i <- i + 1L
        }
        if (mat[i-1,cc]==lowcut) {
            ## fill in more values near the minimum
        }
        if (mat[i-1,cc]==upcut) {
            ## fill in more values near the maximum
        }
        mat
    }

    ## bounds on Cholesky: [0,Inf) for diag, (-Inf,Inf) for off-diag
    ## bounds on sd-corr:  [0,Inf) for diag, (-1.0,1.0) for off-diag
    lower <- pmax(fitted@lower, -1.)
    upper <- 1/(fitted@lower != 0)## = ifelse(fitted@lower==0, Inf, 1.0)
    if (useSc) { # bounds for sigma
        lower <- c(lower,0)
        upper <- c(upper,Inf)
    }
    ## bounds on fixed parameters (TODO: allow user-specified bounds, e.g. for NLMMs)
    lower <- c(lower,rep.int(-Inf, p))
    upper <- c(upper, rep.int(Inf, p))
    npar1 <- if (isLMM(fitted)) nvp else nptot
    ## check that devfun2() computation for the base parameters is (approx.) the
    ##  same as the original devfun() computation
    if(!isTRUE(all.equal(unname(dd(opt[seq(npar1)])), base, tolerance=1e-5))){
        stop("Profiling over both the residual variance and\n",
             "fixed effects is not numerically consistent with\n",
             "profiling over the fixed effects only")}

    ## sequence of variance parameters to profile
    seqnvp <- intersect(seq_len(npar1), which)
    ## sequence of 'all' parameters
    seqpar1 <- seq_len(npar1)
    lowvp <- lower[seqpar1]
    upvp <- upper[seqpar1]
    form <- .zeta ~ foo # pattern for interpSpline formula

    ## as in bootMer.R, define FUN as a
    ##    closure containing the referenced variables
    ##    in its scope to avoid explicit clusterExport statement
    ##    in the PSOCKcluster case
    FUN <- local({
        function(w) {
        if (verbose) cat(if(isLMM(fitted)) "var-cov " else "",
                             "parameter ",w,":\n",sep="")
        wp1 <- w + 1L
        start <- opt[seqpar1][-w]
        pw <- opt[w]
        lowcut <- lower[w]
        upcut <- upper[w]
        zeta <- function(xx,start) {
            ores <- tryCatch(optwrap(optimizer, par=start,
                                     fn=function(x) dd(mkpar(npar1, w, xx, x)),
                                     lower = lowvp[-w],
                                     upper = upvp [-w],
                                     control = control),
                             error=function(e) NULL)
            if (is.null(ores)) {
                devdiff <- NA
                pars <- NA
            } else {
                devdiff <- ores$fval - base
                pars <- ores$par
            }
            if (is.na(devdiff)) {
                warning("NAs detected in profiling")
            } else {
                if(verbose && devdiff < 0)
                    cat("old deviance ",base,",\n",
                        "new deviance ",ores$fval,",\n",
                        "new params ",
                        paste(mkpar(npar1,w,xx,ores$par),
                              collapse=","),"\n")
                if (devdiff < (-devtol))
                    stop("profiling detected new, lower deviance")
                if(devdiff < 0)
                    warning(gettextf("slightly lower deviances (diff=%g) detected",
                                     devdiff), domain=NA)
            }
            devdiff <- max(0,devdiff)
            zz <- sign(xx - pw) * sqrt(devdiff)
            r <- c(zz, mkpar(npar1, w, xx, pars))
            if (isLMM(fitted)) c(r, pp$beta(1)) else r
        }## {zeta}

        ## intermediate storage for pos. and neg. increments
        pres <- nres <- res
        ## assign one row, determined by inc. sign, from a small shift
        ## FIXME:: do something if pw==0 ???
        shiftpar <- if (pw==0) 1e-3 else pw*1.01
        ## Since both the pos- and neg-increment matrices are already
        ## filled with the opt. par. results, this sets the first
        ## two rows of the positive-increment matrix
        ## to (opt. par, shiftpar) and the first two rows of
        ## the negative-increment matrix to (shiftpar, opt. par),
        ## which sets up two points going in the right direction
        ## for each matrix (since the profiling algorithm uses increments
        ## between rows to choose the next parameter increment)
        nres[1, ] <- pres[2, ] <- zeta(shiftpar, start=opt[seqpar1][-w])
        ## fill in the rest of the arrays and collapse them
        upperf <- fillmat(pres, lowcut, upcut, zeta, wp1)
        lowerf <-
            if (pw > lowcut)
                fillmat(nres, lowcut, upcut, zeta, wp1)
            else ## don't bother to fill in 'nres' matrix
                nres
        ## this will throw away the extra 'opt. par' and 'shiftpar'
        ## rows introduced above:
        bres <- as.data.frame(unique(rbind2(upperf,lowerf)))
        pname <- names(opt)[w]
        bres$.par <- pname
        bres <- bres[order(bres[, wp1]), ]

        ## FIXME: test for bad things here??
        form[[3]] <- as.name(pname)
        forspl <- NULL # (in case of error)
        ## bakspl
	bakspl <-
	    tryCatch(backSpline(
		forspl <- interpSpline(form, bres, na.action=na.omit)),
		     error=function(e)e)
        if (inherits(bakspl, "error"))
            warning("non-monotonic profile for ",pname)
        ## return:
        namedList(bres,bakspl,forspl) # namedList() -> lmerControl.R

    }}) ## FUN()

    ## copied from bootMer: DRY!
    L <- if (do_parallel) {
        if (have_mc) {
            parallel::mclapply(seqnvp, FUN, mc.cores = ncpus)
        } else if (have_snow) {
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                ## explicit export of the lme4 namespace since most FUNs will probably
                ## use some of them
                parallel::clusterExport(cl, varlist=getNamespaceExports("lme4"))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl)
                pres <- parallel::parLapply(cl, seqnvp, FUN)
                parallel::stopCluster(cl)
                pres
            } else parallel::parLapply(cl, seqnvp, FUN)
        }
    } else lapply(seqnvp, FUN)
    nn <- names(opt[seqnvp])
    ans <-    setNames(lapply(L,`[[`,"bres"),nn)
    bakspl <- setNames(lapply(L,`[[`,"bakspl"),nn)
    forspl <- setNames(lapply(L,`[[`,"forspl"),nn)

    ## profile fixed effects separately (for LMMs)
    if (isLMM(fitted)) {
        offset.orig <- fitted@resp$offset
        fp <- seq_len(p)
        fp <- fp[(fp+nvp) %in% which]
        ## FIXME: parallelize this too ...
        for (j in fp) {
            if (verbose) cat("fixed-effect parameter ",j,":\n",sep="")
            pres <-            # intermediate results for pos. incr.
                nres <- res    # and negative increments
            est <- opt[nvp + j]
            std <- stderr[j]
            Xw <- X.orig[, j, drop=TRUE]
            Xdrop <- .modelMatrixDrop(X.orig, j)
            pp1 <- new(class(pp),
                       X = Xdrop,
                       Zt = pp$Zt,
                       Lambdat = pp$Lambdat,
                       Lind = pp$Lind,
                       theta = pp$theta,
                       n = nrow(Xdrop))
### FIXME Change this to use the deep copy and setWeights, setOffset, etc.
            rr <- new(Class=class(fitted@resp), y=fitted@resp$y)
            rr$setWeights(fitted@resp$weights)
            fe.zeta <- function(fw, start) {
                ## (start parameter ignored)
                rr$setOffset(Xw * fw + offset.orig)
		rho <- list2env(list(pp=pp1, resp=rr), parent = parent.frame())
		ores <- optwrap(optimizer, par = thopt, fn = mkdevfun(rho, 0L),
				lower = fitted@lower)
                ## this optimization is done on the ORIGINAL
                ##   theta scale (i.e. not the sigma/corr scale)
                ## upper=Inf for all cases
                ## lower = pmax(fitted@lower, -1.0),
                ## upper = 1/(fitted@lower != 0))## = ifelse(fitted@lower==0, Inf, 1.0)
                fv <- ores$fval
                sig <- sqrt((rr$wrss() + pp1$sqrL(1))/n)
                c(sign(fw - est) * sqrt(fv - base),
                  Cv_to_Sv(ores$par, lengths(fitted@cnms), s=sig),
                  ## ores$par * sig, sig,
                  mkpar(p, j, fw, pp1$beta(1)))
            }
            nres[1, ] <- pres[2, ] <- fe.zeta(est + delta * std)
            poff <- nvp + 1L + j
            ## Workaround R bug [rbind2() is S4 generic; cannot catch warnings in its arg]
            ## see lme4 GH issue #304
            upperf <- fillmat(pres, -Inf, Inf, fe.zeta, poff)
            lowerf <- fillmat(nres, -Inf, Inf, fe.zeta, poff)
            bres <- as.data.frame(unique(rbind2(upperf, lowerf)))
            bres$.par <- n.j <- names(fe.orig)[j]
            ans[[n.j]] <- bres[order(bres[, poff]), ]
            form[[3]] <- as.name(n.j)
            bakspl[[n.j]] <-
                tryCatch(backSpline(forspl[[n.j]] <- interpSpline(form, bres)),
                         error=function(e)e)
            if (inherits(bakspl[[n.j]], "error"))
                warning("non-monotonic profile for ", n.j)
        } ## for(j in 1..p)
    } ## if isLMM

    ans <- do.call(rbind, ans)
    row.names(ans) <- NULL ## TODO: rbind(*, make.row.names=FALSE)
    ans$.par <- factor(ans$.par)
    structure(ans,
	      forward = forspl,
	      backward = bakspl,
	      lower = lower[seqnvp],
	      upper = upper[seqnvp],
	      class = c("thpr", "data.frame"))# 'thpr': see ../man/profile-methods.Rd
} ## profile.merMod

##' Transform 'which' \in {parnames | integer | "beta_" | "theta_"}
##' into integer indices
##' @param which numeric or character vector
##' @param nvp number of variance-covariance parameters
##' @param nptot total number of parameters
##' @param parnames vector of parameter names
##' @param verbose print messages?
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
##' tn <- names(getME(fm1,"theta"))
##' nn <- c(tn,names(fixef(fm1)))
##' get.which("theta_",length(tn),length(nn),nn, verbose=TRUE)
##'
get.which <- function(which, nvp, nptot, parnames, verbose=FALSE) {
    if (is.null(which))
        seq_len(nptot)
    else if (is.character(which)) {
        wi <- integer(); wh <- which
        if(any(j <- wh == "theta_")) {
            wi <- seq_len(nvp); wh <- wh[!j] }
        if(any(j <- wh == "beta_") && nptot > nvp) {
            wi <- c(wi, seq(nvp+1, nptot)); wh <- wh[!j] }
        if(any(j <- parnames %in% wh)) { ## which containing param.names
            wi <- sort(unique(c(wi, seq_len(nptot)[j])))
        }
        if(verbose) message(gettextf("From original which = %s: new which <- %s",
                                     deparse(which, nlines=1), deparse(wi, nlines=1)),
                            domain=NA)
        if(length(wi) == 0)
            warning(gettextf("Nothing selected by 'which=%s'", deparse(which)),
                    domain=NA)
        wi
    } else #  stopifnot( .. numeric ..)
        which
}

## This is a hack.  The preferred approach is to write a
## subset method for the ddenseModelMatrix and dsparseModelMatrix
## classes
.modelMatrixDrop <- function(mm, w) {
    if (isS4(mm)) {
	nX <- slotNames(X <- mm[, -w, drop = FALSE])
	do.call(new,
		c(list(Class = class(mm),
		       assign = attr(mm,"assign")[-w],
		       contrasts = NULL
		       ## FIXME: where did the contrasts information go??
		       ##      mm@contrasts
		       ),
		  lapply(structure(nX, .Names=nX),
			 function(nm) slot(X, nm))))
    } else {
	structure(mm[, -w, drop=FALSE],
		  assign = attr(mm, "assign")[-w])
    }
}

## The deviance is profiled with respect to the fixed-effects
## parameters but not with respect to sigma. The other parameters
## are on the standard deviation scale, not the theta scale.
##
## @title Return a function for evaluation of the deviance.
## @param fm a fitted model of class merMod
## @return a function for evaluating the deviance in the extended
##     parameterization.  This is profiled with respect to the
##     variance-covariance parameters (fixed-effects done separately).
devfun2 <- function(fm, useSc, signames)
{
    ## FIXME: have to distinguish between
    ## 'useSc' (GLMM: report profiled scale parameter) and
    ## 'useSc' (NLMM/LMM: scale theta by sigma)
    ## GLMMuseSc <- fm@devcomp$dims["useSc"]
    stopifnot(is(fm, "merMod"))
    fm <- refitML(fm)
    basedev <- -2*c(logLik(fm))  ## no longer deviance()
    vlist <- lengths(fm@cnms)
    sig <- sigma(fm)  ## only if useSc=TRUE?
    stdErr <- unname(coef(summary(fm))[,2])
    pp <- fm@pp$copy()
    ## opt <- c(pp$theta*sig, sig)
    if (useSc) {
        opt <- Cv_to_Sv(pp$theta, n=vlist, s=sig)
    } else {
        opt <- Cv_to_Sv(pp$theta, n=vlist)
    }
    names(opt) <- profnames(fm, signames)
    opt <- c(opt, fixef(fm))
    resp <- fm@resp$copy()
    np <- length(pp$theta)
    nf <- length(fixef(fm))
    if (!isGLMM(fm)) np <- np + 1L
    n <- nrow(pp$V)                   # use V, not X so it works with nlmer
    if (isLMM(fm)) {
        ans <- function(pars)
        {
            stopifnot(is.numeric(pars), length(pars) == np)
            ## Assumption:  all parameters, including the residual SD on SD-scale
            sigma <- pars[np]
            ## .Call(lmer_Deviance, pp$ptr(), resp$ptr(), pars[-np]/sigma)
            ## convert from sdcor vector back to 'unscaled theta'
            thpars <- Sv_to_Cv(pars,n=vlist,s=sigma)
            .Call(lmer_Deviance, pp$ptr(), resp$ptr(), thpars)
            sigsq <- sigma^2
            pp$ldL2() - ldW + (resp$wrss() + pp$sqrL(1))/sigsq + n * log(2 * pi * sigsq)
        }
        ldW <- sum(log(environment(ans)$resp$weights))
        assign("ldW", ldW, envir = environment(ans))
    } else {
        d0 <- update(fm,devFunOnly=TRUE)
        ## from glmer:
        ## rho <- new.env(parent=parent.env(environment()))
        ## rho$pp <- do.call(merPredD$new, c(reTrms[c("Zt","theta","Lambdat","Lind")], n=nrow(X), list(X=X)))
        ## rho$resp <- mkRespMod(fr, if(REML) p else 0L)
        ans <- function(pars)
        {
            stopifnot(is.numeric(pars), length(pars) == np+nf)
            ## FIXME: allow useSc (i.e. NLMMs)
            if (!useSc) {
                thpars <- Sv_to_Cv(pars[seq(np)],n=vlist)
            } else {
                thpars <- Sv_to_Cv(pars[seq(np)],n=vlist,s=pars[np])
            }
            fixpars <- pars[-seq(np)]
            d0(c(thpars,fixpars))
        }
    }
    attr(ans, "optimum") <- opt         # w/ names()
    attr(ans, "basedev") <- basedev
    attr(ans, "thopt") <- pp$theta
    attr(ans, "stderr") <- stdErr
    class(ans) <- "devfun"
    ans
}

## extract only the y component from a prediction
predy <- function(sp, vv) {
    if (inherits(sp,"error")) rep(NA,length(vv))
    else predict(sp, vv)$y
}

stripExpr <- function(ll, nms) {
    stopifnot(is.list(ll), is.character(nms))
    fLevs <- as.expression(nms) # variable names; now replacing log(sigma[.]) etc:
    fLevs[nms ==  ".sigma"] <- expression(sigma)
    fLevs[nms == ".lsigma"] <- expression(log(sigma))
    fLevs[nms == ".sigmasq"] <- expression(sigma^2)
    sigNms  <- grep("^\\.sig[0-9]+", nms)
    lsigNms <- grep("^\\.lsig[0-9]+", nms)
    sig2Nms <- grep("^\\.sigsq[0-9]+", nms)
    ## the <n> in ".sig0<n>" and then in ".lsig0<n>":
    sigsub  <- as.integer(substring(nms[ sigNms], 5))
    lsigsub <- as.integer(substring(nms[lsigNms], 6))
    sig2sub <- as.integer(substring(nms[sig2Nms], 7))
    fLevs[ sigNms] <- lapply( sigsub, function(i) bquote(    sigma[.(i)]))
    fLevs[lsigNms] <- lapply(lsigsub, function(i) bquote(log(sigma[.(i)])))
    fLevs[sig2Nms] <- lapply(sig2sub, function(i) bquote(   {sigma[.(i)]}^2))
    ## result of using { .. }^2  is easier to understand    ==          ==
    levsExpr <- substitute(strip.custom(factor.levels=foo), list(foo=fLevs))
    llNms <- names(ll)
    snames <- c("strip", "strip.left")
    if (all(!(snames %in% llNms))) {
        ll$strip <- levsExpr
    } else {
        lapply(snames, function(nm) {
            if (nm %in% llNms) {
                vv <- ll[[nm]]
                if (is.logical(vv) && vv) ll[[nm]] <<- levsExpr
            }
        })
    }
    ll
}

panel.thpr <- function(x, y, spl, absVal, ...)
{
    panel.grid(h = -1, v = -1)
    myspl <- spl[[panel.number()]]
    lsegments(x, y, x, 0, ...)
    if (absVal) {
        y[y == 0] <- NA
        lsegments(x, y, rev(x), y)
    } else {
        panel.abline(h = 0, ...)
    }
    if (!is(myspl,"spline")) {
        ## 'difficult' data
        if (absVal) myspl$y <- abs(myspl$y)
        panel.lines (myspl$x, myspl$y)
        panel.points(myspl$x, myspl$y, pch="+")
        warning(gettextf("bad profile for variable %d: using linear interpolation",
                         panel.number()), domain=NA)
    } else {
        lims <- current.panel.limits()$xlim
        krange <- range(myspl$knots)
        pr <- predict(myspl,
                      seq(max(lims[1], krange[1]),
                          min(lims[2], krange[2]), len = 101))
        if (absVal) pr$y <- abs(pr$y)
        panel.lines(pr$x, pr$y)
    }
}

## A lattice-based plot method for profile objects
##' @importFrom lattice xyplot
##' @S3method xyplot thpr
xyplot.thpr <-
    function (x, data = NULL,
              levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), df = 1)),
              conf = c(50, 80, 90, 95, 99)/100,
              absVal = FALSE, scales = NULL,
              which = 1:nptot, ...)
{
    if(any(!is.finite(conf) | conf <= 0 | conf >= 1))
        stop("values of 'conf' must be strictly between 0 and 1")
    stopifnot(1 <= (nptot <- length(nms <- levels(x[[".par"]]))))
    ## FIXME: is this sufficiently reliable?
    ## (include "sigma" in 'theta' parameters)
    nvp <- length(grep("^(\\.sig[0-9]+|.sigma|sd_|cor_)", nms))
    which <- get.which(which, nvp, nptot, nms)
    levels <- sort(levels[is.finite(levels) & levels > 0])
    spl  <- attr(x, "forward") [which]
    bspl <- attr(x, "backward")[which]
    ## for parameters for which spline couldn't be computed,
    ## replace the 'spl' element with the raw profile data
    if(any(badSpl <- vapply(spl, is.null, NA))) {
	spl[badSpl] <- lapply(which(badSpl), function(i) {
	    n <- names(badSpl)[i]
	    r <- x[x[[".par"]] == n, ]
	    data.frame(y = r[[".zeta"]], x = r[[n]])
	})
	bspl[badSpl] <- lapply(spl[badSpl], function(d) data.frame(x=d$y,y=d$x))
	## FIXME: more efficient, not yet ok ?
	## ibad <- which(badSpl)
	## spl[ibad] <- lapply(names(ibad), function(n) {
	##     r <- x[x[[".par"]]==n,]
	##     data.frame(y = r[[".zeta"]], x = r[[n]])
	## })
	## bspl[ibad] <- lapply(spl[ibad], function(d) data.frame(x=d$y,y=d$x))
    }
    zeta <- c(-rev(levels), 0, levels)
    mypred <- function(bs, zeta) { ## use linear approximation if backspline doesn't work
	if (inherits(bs,"spline"))
	    predy(bs, zeta)
	else if(is.numeric(x <- bs$x) && is.numeric(y <- bs$y) && length(x) == length(y))
	    approx(x, y, xout = zeta)$y
	else
	    rep_len(NA, length(zeta))
    }
    fr <- data.frame(zeta = rep.int(zeta, length(spl)),
                     pval = unlist(lapply(bspl, mypred, zeta)),
                     pnm = gl(length(spl), length(zeta), labels = names(spl)))
    if (length(ind <- which(is.na(fr$pval)))) {
        fr[ind, "zeta"] <- 0
        for (i in ind)
### FIXME: Should check which bound has been violated, although it
### will almost always be the minimum.
            if (inherits(curspl <-  spl[[fr[i, "pnm"] ]], "spline")) {
                fr[i, "pval"] <- min(curspl$knots)
            }
    }
    ylab <- if (absVal) {
        fr$zeta <- abs(fr$zeta)
        expression("|" * zeta * "|")
    } else
        expression(zeta)
    intscales <- list(x = list(relation = 'free'))
    ## FIXME: is there something less clunky we can do here
    ##   that allows for all possible user inputs
    ##   (may want to (1) override x$relation (2) add elements to $x
    ##    (3) add elements to scales)
    if (!is.null(scales)) {
        if (!is.null(scales[["x"]])) {
            if (!is.null(scales[["x"]]$relation)) {
                intscales[["x"]]$relation <- scales[["x"]]$relation
                scales[["x"]]$relation <- NULL
            }
            intscales[["x"]] <- c(intscales[["x"]],scales[["x"]])
            scales[["x"]] <- NULL
        }
        intscales <- c(intscales,scales)
    }
    ll <- c(list(...),
            list(x = zeta ~ pval | pnm, data=fr,
                 scales = intscales,
                 ylab = ylab, xlab = NULL, panel=panel.thpr,
                 spl = spl, absVal = absVal))

    do.call(xyplot, stripExpr(ll, names(spl)))
}

## copy of stats:::format.perc (not exported, and ":::" being forbidden nowadays):
format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE,
                 scientific = FALSE, digits = digits),
          "%")
}

##' confint() method for  our profile() results 'thpr'
##' @importFrom stats confint
confint.thpr <- function(object, parm, level = 0.95, zeta,
                         ## tolerance for non-monotonic profiles
                         ## (raw values, not splines)
                         non.mono.tol=1e-2,
                         ...)
{
    bak <- attr(object, "backward")
    ## fallback strategy for old profiles that don't have a lower/upper
    ##  attribute saved ...
    if (is.null(lower <- attr(object,"lower")))
        lower <- rep(NA,length(parm))
    if (is.null(upper <- attr(object,"upper")))
        upper <- rep(NA,length(parm))
    ## FIXME: work a little harder to add -Inf/Inf for fixed effect
    ##  parameters?  (Should only matter for really messed-up profiles)
    bnms <- names(bak)
    if (missing(parm)) parm <- bnms
    else if (is.numeric(parm)) parm <- bnms[parm]
    parm <- intersect(as.character(parm), bnms)
    cn <- NULL
    if (missing(zeta)) {
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        zeta <- qnorm(a)
        cn <- format.perc(a, 3)
    }
    ci <- matrix(NA,nrow=length(parm),ncol=2,
                 dimnames=list(parm,cn))
    for (i in seq_along(parm)) {
        ## would like to build this machinery into predy, but
        ## predy is used in many places and it's much harder to
        ## tell in general whether an NA indicates a lower or an
        ## upper bound ...
        badprof <- FALSE
        p <- rep(NA,2)
	if (!inherits(b <- bak[[parm[i]]], "error")) {
            p <- predy(b, zeta)
        } else {
            obj1 <- object[object$.par==parm[[i]],c(parm[[i]],".zeta")]
            if (all(is.na(obj1[,2]))) {
                badprof <- TRUE
                warning("bad profile for ",parm[i])
            } else if (min(diff(obj1[,2])<(-non.mono.tol),na.rm=TRUE)) {
                badprof <- TRUE
                warning("non-monotonic profile for ",parm[i])
            } else {
                warning("bad spline fit for ",parm[i],": falling back to linear interpolation")
                p <- approxfun(obj1[,2],obj1[,1])(zeta)
            }
        }
        if (!badprof) {
            if (is.na(p[1])) p[1] <- lower[i]
            if (is.na(p[2])) p[2] <- upper[i]
        }
        ci[i,] <- p
    }
    ci
}

## FIXME: make bootMer more robust; make profiling more robust;
## more warnings; documentation ...

##' Compute confidence intervals on the parameters of an lme4 fit
##' @param object a fitted [ng]lmer model
##' @param parm parameters (specified by integer position)
##' @param level confidence level
##' @param method for computing confidence intervals
##' @param zeta likelihood cutoff
##' (if not specified, computed from \code{level}: "profile" only)
##' @param nsim number of simulations for parametric bootstrap intervals
##' @param boot.type bootstrap confidence interval type
##' @param quiet (logical) suppress messages about computationally intensive profiling?
##' @param oldNames (logical) use old-style names for \code{method="profile"}? (See \code{signames} argument to \code{\link{profile}}
##' @param \dots additional parameters to be passed to  \code{\link{profile.merMod}} or \code{\link{bootMer}}
##' @return a numeric table of confidence intervals
confint.merMod <- function(object, parm, level = 0.95,
			   method = c("profile","Wald","boot"),
			   zeta, nsim=500, boot.type = c("perc","basic","norm"),
                           FUN = NULL, quiet=FALSE, oldNames=TRUE, ...)
{
    method <- match.arg(method)
    boot.type <- match.arg(boot.type)
    ## 'parm' corresponds to 'which' in other contexts
    if (method=="boot" && !is.null(FUN)) {
        ## custom boot function, don't expand parameter names
    } else {
        vn <- profnames(object,oldNames)
        an <- c(vn,names(fixef(object)))
        if (missing(parm)) {
            parm <- seq(length(an))
        } else {
            parm <- get.which(parm, nvp=length(vn), nptot=length(an),
                              parnames=an)
        }
        if (!quiet && method %in% c("profile","boot")) {
            mtype <- switch(method, profile="profile", boot="bootstrap")
            message("Computing ",mtype," confidence intervals ...")
            flush.console()
        }
    }
    switch(method,
	   "profile" =
           {
               pp <- profile(object, which=parm, signames=oldNames, ...)
               confint(pp, level=level, zeta=zeta)
           },
	   "Wald" =
           {
               a <- (1 - level)/2
               a <- c(a, 1 - a)
               ci.vcov <- array(NA,dim = c(length(vn), 2L),
                                dimnames = list(vn, format.perc(a,3)))
               ## copied with small changes from confint.default
               cf <- fixef(object)   ## coef() -> fixef()
               pnames <- names(cf)
               ## N.B. can't use sqrt(...)[parm] (diag() loses names)
               ses <- sqrt(diag(vcov(object)))
               ci.fixed <- array(cf + ses %o% qnorm(a),
                                 dim = c(length(pnames), 2L),
                                 dimnames = list(pnames, format.perc(a, 3)))
               vnames <- tnames(object)
               ci.all <- rbind(ci.vcov,ci.fixed)
               ci.all[parm,,drop=FALSE]
           },
	   "boot" =
           {
               bootFun <- function(x) {
		   th <- x@theta
		   nvec <- lengths(x@cnms)
                   scaleTh <- (isLMM(x) || isNLMM(x))
                   useSc <- as.logical(x@devcomp$dims[["useSc"]])
		   ## FIXME: still ugly.  Best cleanup via Cv_to_Sv ...
		   ss <- if (scaleTh) {	 ## scale variances by sigma and include it
		       setNames(Cv_to_Sv(th,n=nvec,s=sigma(x)), vn)
		   } else if (useSc) { ## don't scale variances but do include sigma
		       setNames(c(Cv_to_Sv(th,n=nvec),sigma(x)), vn)
		   } else {  ## no scaling, no sigma
		       setNames(Cv_to_Sv(th,n=nvec), vn)
		   }
                   c(ss, fixef(x))
               }
               if (is.null(FUN)) FUN <- bootFun
               bb <- bootMer(object, FUN=FUN, nsim=nsim,...)
               bci <- lapply(seq_along(bb$t0),
                             boot.out=bb,
                             boot::boot.ci, type=boot.type, conf=level)
               cpos <- grep(boot.type, names(bci[[1]]))
               ## get _last_ two columns
               ccol <- ncol(bci[[1]][[cpos]])+(-1:0)
               citab <- t(sapply(bci,function(x) x[[cpos]][ccol]))
               a <- (1 - level)/2
               a <- c(a, 1 - a)
               dimnames(citab) <- list(names(bb[["t0"]]), format.perc(a, 3))
               if (missing(parm)) {
                   ## only happens if we have custom boot method
                   parm <- rownames(citab)
               }
               citab[parm, , drop=FALSE]
           },
           ## should never get here ...
           stop("unknown confidence interval method"))
}

##' Convert x-cosine and y-cosine to average and difference.
##'
##' Convert the x-cosine and the y-cosine to an average and difference
##' ensuring that the difference is positive by flipping signs if
##' necessary
##' @param xc x-cosine
##' @param yc y-cosine
ad <- function(xc, yc)
{
    a <- (xc + yc)/2
    d <- (xc - yc)
    cbind(sign(d)* a, abs(d))
}

##' convert d versus a (as an xyVector) and level to a matrix of taui and tauj
##' @param xy an xyVector
##' @param lev the level of the contour
tauij <- function(xy, lev) lev * cos(xy$x + outer(xy$y/2, c(-1, 1)))

##' @title safe arc-cosine
##' @param x numeric vector argument
##' @return acos(x) being careful of boundary conditions
sacos <- function(x) acos(pmax.int(-0.999, pmin.int(0.999, x)))

##' Generate a contour
##'
##' @title Generate a contour
##' @param sij the arc-cosines of i on j
##' @param sji the arc-cosines of j on i
##' @param levels numeric vector of levels at which to interpolate
##' @param nseg number of segments in the interpolated contour
##' @return a list with components
##' \item{tki}{the tau-scale predictions of i on j at the contour levels}
##' \item{tkj}{the tau-scale predictions of j on i at the contour levels}
##' \item{pts}{an array of dimension (length(levels), nseg, 2) containing the points on the contours}
cont <- function(sij, sji, levels, nseg = 101)
{
    ada <- array(0, c(length(levels), 2, 4))
    ada[, , 1] <- ad(0, sacos(predy(sij,  levels)/levels))
    ada[, , 2] <- ad(sacos(predy(sji, levels)/levels), 0)
    ada[, , 3] <- ad(pi, sacos(predy(sij, -levels)/levels))
    ada[, , 4] <- ad(sacos(predy(sji, -levels)/levels), pi)
    pts <- array(0, c(length(levels), nseg + 1, 2))
    for (i in seq_along(levels))
        pts[i, ,] <- tauij(predict(periodicSpline(ada[i, 1, ], ada[i, 2, ]),
                                   nseg = nseg), levels[i])
    levs <- c(-rev(levels), 0, levels)
    list(tki = predict(sij, levs), tkj = predict(sji, levs), pts = pts)
}

## copied from lattice:::chooseFace
chooseFace <- function (fontface = NULL, font = 1)
{
    if (is.null(fontface))
        font
    else fontface
}


##' Draws profile pairs plots.  Contours are for the marginal
##' two-dimensional regions (i.e. using df = 2).
##'
##' @title Profile pairs plot
##' @param x the result of \code{\link{profile}} (or very similar structure)
##' @param data unused - only for compatibility with generic
##' @param levels the contour levels to be shown; usually derived from \code{conf}
##' @param conf numeric vector of confidence levels to be shown as contours
##' @param ... further arguments passed to \code{\link{splom}}
##' @importFrom grid gpar viewport
##' @importFrom lattice splom
##' @method splom thpr
##' @export
splom.thpr <- function (x, data,
                        levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 2)),
                        conf = c(50, 80, 90, 95, 99)/100, which = 1:nptot,
                        draw.lower = TRUE, draw.upper = TRUE, ...)
{
    stopifnot(1 <= (nptot <- length(nms <- names(attr(x, "forward")))))
    singfit <- FALSE
    for (i in grep("^(\\.sig[0-9]+|sd_)", names(x)))
        singfit <- singfit || any(x[,".zeta"] == 0  &  x[,i] == 0)
    if (singfit) warning("splom is unreliable for singular fits")

    nvp <- length(grep("^(\\.sig[0-9]+|.sigma|sd_|cor_)", nms))
    which <- get.which(which, nvp, nptot, nms)
    if (length(which) == 1)
        stop("can't draw a scatterplot matrix for a single variable")

    mlev <- max(levels)
    spl <- attr(x, "forward")[which]
    frange <- sapply(spl, function(x) range(x$knots))
    bsp <- attr(x, "backward")[which]
    x <- x[x[[".par"]] %in% nms[which],c(".zeta",nms[which],".par")]
    ## brange <- sapply(bsp, function(x) range(x$knots))
    pfr <- do.call(cbind, lapply(bsp, predy, c(-mlev, mlev)))
    pfr[1, ] <- pmax.int(pfr[1, ], frange[1, ], na.rm = TRUE)
    pfr[2, ] <- pmin.int(pfr[2, ], frange[2, ], na.rm = TRUE)
    nms <- names(spl)
    ## Create data frame fr of par. vals in zeta coordinates
    fr <- x[, -1]
    for (nm in nms) fr[[nm]] <- predy(spl[[nm]], na.omit(fr[[nm]]))
    fr1 <- fr[1, nms]
    ## create a list of lists with the names of the parameters
    traces <- lapply(fr1, function(el) lapply(fr1, function(el1) list()))
    for (j in seq_along(nms)[-1]) {
        for (i in seq_len(j - 1)) {
            .par <- NULL  ## suppress R CMD check warning
            fri <- subset(fr, .par == nms[i])
            sij <- interpSpline(fri[ , i], fri[ , j])
            frj <- subset(fr, .par == nms[j])
            sji <- interpSpline(frj[ , j], frj[ , i])
            ll <- cont(sij, sji, levels)
            traces[[j]][[i]] <- list(sij = sij, sji = sji, ll = ll)
        }
    }
    if(draw.lower) ## panel function for lower triangle
    lp <- function(x, y, groups, subscripts, i, j, ...) {
        tr <- traces[[j]][[i]]
        grid::pushViewport(viewport(xscale = c(-1.07, 1.07) * mlev,
                                    yscale = c(-1.07, 1.07) * mlev))
        dd <- sapply(current.panel.limits(), diff)/50
        psij <- predict(tr$sij)
        ll <- tr$ll
        ## now do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(psij$y, psij$x, ...)
        llines(predict(tr$sji), ...)
        with(ll$tki, lsegments(y - dd[1], x, y + dd[1], x, ...))
        with(ll$tkj, lsegments(x, y - dd[2], x, y + dd[2], ...))
        for (k in seq_along(levels)) llines(ll$pts[k, , ], ...)
        grid::popViewport(1)
    }
    if(draw.upper) ## panel function for upper triangle
    up <- function(x, y, groups, subscripts, i, j, ...) {
        ## panels are transposed so reverse i and j
        jj <- i
        ii <- j
        tr <- traces[[jj]][[ii]]
        ll <- tr$ll
        pts <- ll$pts
        ## limits <- current.panel.limits()
        psij <- predict(tr$sij)
        psji <- predict(tr$sji)
        ## do the actual plotting
        panel.grid(h = -1, v = -1)
        llines(predy(bsp[[ii]], psij$x), predy(bsp[[jj]], psij$y), ...)
        llines(predy(bsp[[ii]], psji$y), predy(bsp[[jj]], psji$x), ...)
        for (k in seq_along(levels))
            llines(predy(bsp[[ii]], pts[k, , 2]),
                   predy(bsp[[jj]], pts[k, , 1]), ...)
    }
    dp <- function(x = NULL,            # diagonal panel
                   varname = NULL, limits, at = NULL, lab = NULL,
                   draw = TRUE,

                   varname.col = add.text$col,
                   varname.cex = add.text$cex,
                   varname.lineheight = add.text$lineheight,
                   varname.font = add.text$font,
                   varname.fontfamily = add.text$fontfamily,
                   varname.fontface = add.text$fontface,

                   axis.text.col = axis.text$col,
                   axis.text.alpha = axis.text$alpha,
                   axis.text.cex = axis.text$cex,
                   axis.text.font = axis.text$font,
                   axis.text.fontfamily = axis.text$fontfamily,
                   axis.text.fontface = axis.text$fontface,

                   axis.line.col = axis.line$col,
                   axis.line.alpha = axis.line$alpha,
                   axis.line.lty = axis.line$lty,
                   axis.line.lwd = axis.line$lwd,
                   i, j,
                   ...)
    {
        n.var <- eval.parent(expression(n.var))
        add.text <- trellis.par.get("add.text")
        axis.line <- trellis.par.get("axis.line")
        axis.text <- trellis.par.get("axis.text")
        if (!is.null(varname))
            grid::grid.text(varname,
                            gp =
                                gpar(col = varname.col,
                                     cex = varname.cex,
                                     lineheight = varname.lineheight,
                                     fontface = chooseFace(varname.fontface,
                                     varname.font),
                                     fontfamily = varname.fontfamily))
        if (draw)
        {
            at <- pretty(limits)
            sides <- c("left", "top")
            if (j == 1) sides <- "top"
            if (j == n.var) sides <- "left"
            for (side in sides)
                panel.axis(side = side,
                           at = at,
                           labels = format(at, trim = TRUE),
                           ticks = TRUE,
                           check.overlap = TRUE,
                           half = side == "top" && j > 1,

                           tck = 1, rot = 0,

                           text.col = axis.text.col,
                           text.alpha = axis.text.alpha,
                           text.cex = axis.text.cex,
                           text.font = axis.text.font,
                           text.fontfamily = axis.text.fontfamily,
                           text.fontface = axis.text.fontface,

                           line.col = axis.line.col,
                           line.alpha = axis.line.alpha,
                           line.lty = axis.line.lty,
                           line.lwd = axis.line.lwd)
            lims <- c(-1.07, 1.07) * mlev
            grid::pushViewport(viewport(xscale = lims, yscale = lims))
	    side <- if(j == 1) "right" else "bottom"
	    which.half <- if(j == 1) "lower" else "upper"
            at <- pretty(lims)
            panel.axis(side = side, at = at, labels = format(at, trim = TRUE),
                       ticks = TRUE, half = TRUE, which.half = which.half,
                       tck = 1, rot = 0,

                       text.col = axis.text.col,
                       text.alpha = axis.text.alpha,
                       text.cex = axis.text.cex,
                       text.font = axis.text.font,
                       text.fontfamily = axis.text.fontfamily,
                       text.fontface = axis.text.fontface,

                       line.col = axis.line.col,
                       line.alpha = axis.line.alpha,
                       line.lty = axis.line.lty,
                       line.lwd = axis.line.lwd)
            grid::popViewport(1)
        }
    }

    panel.blank <- function(...) {}
    splom(~ pfr,
          lower.panel = if(draw.lower) lp else panel.blank,
          upper.panel = if(draw.upper) up else panel.blank,
          diag.panel = dp, ...)
}

## return an lmer profile like x with all the .sigNN parameters
## replaced by .lsigNN.  The forward and backward splines for
## these parameters are recalculated.  -> ../man/profile-methods.Rd
logProf <- function (x, base = exp(1), ranef=TRUE,
                     sigIni = if(ranef) "sig" else "sigma")
{
    stopifnot(inherits(x, "thpr"))
    cn <- colnames(x)
    sigP <- paste0("^\\.", sigIni)
    if (length(sigs <- grep(sigP, cn))) {
        repP <- sub("sig", ".lsig", sigIni)
	colnames(x) <- cn <- sub(sigP, repP, cn)
        levels(x[[".par"]]) <- sub(sigP, repP, levels(x[[".par"]]))
        names(attr(x, "backward")) <-
            names(attr(x, "forward")) <-
                sub(sigP, repP, names(attr(x, "forward")))
	for (nm in cn[sigs]) {
            x[[nm]] <- log(x[[nm]], base = base)
	    fr <- x[x[[".par"]] == nm & is.finite(x[[nm]]), TRUE, drop=FALSE]
            form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
            attr(x, "forward")[[nm]] <- isp <- interpSpline(form, fr)
            attr(x, "backward")[[nm]] <- backSpline(isp)
        }
        ## eliminate rows that produced non-finite logs
        x <- x[apply(is.finite(as.matrix(x[, sigs])), 1, all),]
    }
    x
}
## the log() method must have (x, base); no other arguments
log.thpr <- function (x, base = exp(1)) logProf(x, base=base)




##' Create an approximating density from a profile object
##'
##' @title Approximate densities from profiles
##' @param pr a profile object
##' @param npts number of points at which to evaluate the density
##' @param upper upper bound on cumulative for a cutoff
##' @return a data frame
dens <- function(pr, npts=201, upper=0.999) {
    stopifnot(inherits(pr, "thpr"))
    npts <- as.integer(npts)
    stopifnot(inherits(pr, "thpr"), npts > 0,
              is.numeric(upper), 0.5 < upper, upper < 1)
    spl <- attr(pr, "forward")
    bspl <- attr(pr, "backward")
    zeta <- c(qnorm(1-upper), qnorm(upper))
    rng <- lapply(bspl, function(spl)
    {
        rng <- predy(spl, zeta)
        if (is.na(rng[1])) rng[1] <- 0
        if (is.na(rng[2])) { ## try harder to pick an upper bound
            upper <- 1-10^seq(-4,-1,length=21)
            i <- 1
            while (is.na(rng[2]) && i<=length(upper)) {
                rng[2] <- predy(spl,qnorm(upper[i]))
                i <- i + 1
            }
            if (is.na(rng[2])) {
                warning("can't find an upper bound for the profile")
                return(rep(NA,npts))
            }
        }
        seq(rng[1], rng[2], len=npts)
    })
    fr <- data.frame(pval=unlist(rng),
                     pnm=gl(length(rng), npts, labels=names(rng)))
    dd <- list()
    for (nm in names(rng)) {
        zz <- predy(spl[[nm]], rng[[nm]])
        dd[[nm]] <- dnorm(zz) * predict(spl[[nm]], rng[[nm]], deriv=1)$y
    }
    fr$density <- unlist(dd)
    fr
}

##' Densityplot method for a mixed-effects model profile
##
##' @title densityplot from a mixed-effects profile
##' @param x a mixed-effects profile
##' @param data not used - for compatibility with generic
##' @param ... optional arguments to \code{\link[lattice]{densityplot}()}
##' from package \pkg{lattice}.
##' @return a density plot
##' @examples ## see   example("profile.merMod")
##' @importFrom lattice densityplot
##' @method densityplot thpr
##' @export
densityplot.thpr <- function(x, data, ...) {
    ll <- c(list(...),
            list(x=density ~ pval|pnm,
                 data=dens(x),
                 type=c("l","g"),
                 scales=list(relation="free"),
                 xlab=NULL))
    do.call(xyplot, stripExpr(ll, names(attr(x, "forward"))))
}

##' Transform a mixed-effects profile to the variance scale
varianceProf <- function(x, ranef=TRUE) {
    ## "parallel" to logProf()
    stopifnot(inherits(x, "thpr"))
    cn <- colnames(x)
    if(length(sigs <- grep(paste0("^\\.", if(ranef)"sig" else "sigma"), cn))) {
        ## s/sigma/sigmasq/ ;  s/sig01/sig01sq/	 etc
        sigP <- paste0("^(\\.sig", if(ranef) "(ma)?" else "ma", ")")
	repP <- "\\1sq"
	colnames(x) <- cn <- sub(sigP, repP, cn)
	levels(x[[".par"]]) <- sub(sigP, repP, levels(x[[".par"]]))
	names(attr(x, "backward")) <-
	    names(attr(x, "forward")) <-
		sub(sigP, repP, names(attr(x, "forward")))
	for (nm in cn[sigs]) {
	    x[[nm]] <- x[[nm]]^2
	    ## select rows (and currently drop extra attributes)
	    fr <- x[x[[".par"]] == nm, TRUE, drop=FALSE]
	    form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
	    attr(x, "forward")[[nm]] <- isp <- interpSpline(form, fr)
	    attr(x, "backward")[[nm]] <- backSpline(isp)
        }
    }
    x
}

## convert profile to data frame, adding a .focal parameter to simplify lattice/ggplot plotting
##' @method as.data.frame thpr
##' @param x the result of \code{\link{profile}} (or very similar structure)
##' @export
##' @rdname profile-methods
as.data.frame.thpr <- function(x,...) {
    class(x) <- "data.frame"
    m <- as.matrix(x[,seq(ncol(x))-1]) ## omit .par
    x.p <- x[[".par"]]
    x[[".focal"]] <- m[cbind(seq(nrow(x)),match(x.p,names(x)))]
    x[[".par"]] <- factor(x.p, levels=unique(as.character(x.p))) ## restore order
    x
}


