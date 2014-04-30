## --> ../man/profile-methods.Rd

##' @importFrom splines backSpline interpSpline periodicSpline
##' @importFrom stats profile
##' @method profile merMod
##' @export
profile.merMod <- function(fitted, which=1:nptot, alphamax = 0.01,
			   maxpts = 100, delta = cutoff/8,
                           verbose=0, devtol=1e-9,
                           maxmult = 10,
                           startmethod = "prev",
                           optimizer = "bobyqa",
                           signames = TRUE, ...)
{
  ## FIXME: allow choice of nextstep/nextstart algorithm?
  ## FIXME: by default, get optimizer from within fitted object
  ## FIXME: allow selection of individual variables to profile by name?
  ## FIXME: allow for failure of bounds (non-pos-definite correlation matrices) when >1 cor parameter
  ## FIXME: generalize to GLMMs
  ## (use different devfun;
  ##  be careful with scale parameter;
  ##  profile all parameters at once rather than RE first and then fixed)

    useSc <- isLMM(fitted) || isNLMM(fitted)
    dd <- devfun2(fitted,useSc,signames)
    ## FIXME: figure out to what do here ...
    if (isGLMM(fitted) && fitted@devcomp$dims["useSc"])
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

    which <- get.which(which,nvp,nptot,names(opt),verbose)

    ans <- lapply(opt[which], function(el) NULL)
    bakspl <- forspl <- ans

    res <- c(.zeta = 0, opt)
    res <- matrix(res, nrow = maxpts, ncol = length(res),
                  dimnames = list(NULL, names(res)), byrow = TRUE)
    ## FIXME: why is cutoff based on nptot (i.e. boundary of simultaneous LRT conf region for nptot values)
    ##  when we are computing (at most) 2-parameter profiles here?

    cutoff <- sqrt(qchisq(1 - alphamax, nptot))

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
            if (r>1) {
                if (abs(step) > (maxstep <- abs(maxmult*num))) {
                    maxstep <- sign(step)*maxstep
                    if (verbose) cat(sprintf("capped step at %1.2f (multiplier=%1.2f > %1.2f)\n",
                                             maxstep,abs(step/num),maxmult))
                    step <- maxstep
                }
            }
        }
        min(upper, max(lower, pvals[2] + sign(num) * step))
      }

    nextstart <- function(mat, pind, r, method="global") {
      ## FIXME: indexing may need to be checked (e.g. for fixed-effect parameters)
      switch(method,
             global=opt[seqpar1][-pind],  ## address opt, no zeta column
             prev=mat[r,1+seqpar1][-pind],
             extrap=stop("stub")) ## do something with mat[r-(1:0),1+seqnvp])[-pind] ...
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
### FIXME:  add code to evaluate more rows near the minimum if that
###        constraint was active.
    fillmat <- function(mat, lowcut, upcut, zetafun, cc) {
        nr <- nrow(mat)
        i <- 2L
        while (i < nr && mat[i, cc] > lowcut && mat[i,cc] < upcut &&
               (is.na(curzeta <- abs(mat[i, ".zeta"])) || curzeta <= cutoff)) {
            np <- nextpar(mat, cc, i, delta, lowcut, upcut)
            ns <- nextstart(mat, cc-1, i, startmethod)
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
    lower <- pmax(fitted@lower,-1.0)
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
    seqnvp <- intersect(seq_len(npar1),which)
    ## sequence of 'all' parameters
    seqpar1 <- seq_len(npar1)
    lowvp <- lower[seqpar1]
    upvp <- upper[seqpar1]
    form <- .zeta ~ foo           # pattern for interpSpline formula

    for (w in seqnvp) {
       if (verbose) cat(if(isLMM(fitted)) "var-cov " else "", "parameter ",w,":\n",sep="")
       wp1 <- w + 1L
       start <- opt[seqpar1][-w]
       pw <- opt[w]
       lowcut <- lower[w]
       upcut <- upper[w]
       zeta <- function(xx,start) {
	   ores <- tryCatch(optwrap(optimizer, par=start,
				    fn=function(x) dd(mkpar(npar1, w, xx, x)),
				    lower = lowvp[-w],
				    upper = upvp [-w]), error=function(e)NULL)
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

### FIXME: The starting values for the conditional optimization should
### be determined from recent starting values, not always the global
### optimum values.

### Can do this most easily by taking the change in the other parameter values at
### the two last points and extrapolating.


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
       if (pw>lowcut) {
           lowerf <- fillmat(nres, lowcut, upcut, zeta, wp1)
       } else {
           ## don't bother to fill in 'nres' matrix
           lowerf <- nres
       }
       ## this will throw away the extra 'opt. par' and 'shiftpar'
       ## rows introduced above:
       bres <- as.data.frame(unique(rbind2(upperf,lowerf)))
       pname <- names(opt)[w]
       bres$.par <- pname
       ans[[pname]] <- bres[order(bres[, wp1]), ]
       form[[3]] <- as.name(pname)

       ## FIXME: test for bad things here??
       bakspl[[pname]] <- tryCatch(backSpline(forspl[[pname]] <-
                                              interpSpline(form, bres,na.action=na.omit)),
                                   error=function(e)e)
       if (inherits(bakspl[[pname]],"error")) {
           warning("non-monotonic profile")
     }

    } ## for(w in ..)

    ## profile fixed effects separately (for LMMs)
    if (isLMM(fitted)) {
        offset.orig <- fitted@resp$offset
        fp <- seq_len(p)
        fp <- fp[(fp+nvp) %in% which]
        for (j in fp) {
            if (verbose) cat("fixed-effect parameter ",j,":\n",sep="")
            pres <-            # intermediate results for pos. incr.
                nres <- res    # and negative increments
            est <- opt[nvp + j]
            std <- stderr[j]
            Xw <-X.orig[, j, drop=TRUE]
            Xdrop <- .modelMatrixDrop(X.orig, j)
            pp1 <- do.call("new", list(Class = class(pp),
                                     X = Xdrop,
                                     Zt = pp$Zt,
                                     Lambdat = pp$Lambdat,
                                     Lind = pp$Lind,
                                     theta = pp$theta,
                                     n = nrow(Xdrop))
                           )
### FIXME Change this to use the deep copy and setWeights, setOffset, etc.
            rr <- new(Class=class(fitted@resp), y=fitted@resp$y)
            rr$setWeights(fitted@resp$weights)
            fe.zeta <- function(fw, start) {
                ## (start parameter ignored)
                rr$setOffset(Xw * fw + offset.orig)
                rho <- as.environment(list(pp=pp1, resp=rr))
                parent.env(rho) <- parent.frame()
                ores <- optwrap(optimizer,
                                par=thopt, fn=mkdevfun(rho, 0L),
                                lower = fitted@lower)
                ## ?? this optimization is done on the ORIGINAL
                ## theta scale (i.e. not the sigma/corr scale ??
                ## upper=Inf for all cases
                ## lower = pmax(fitted@lower, -1.0),
                ## upper = 1/(fitted@lower != 0))## = ifelse(fitted@lower==0, Inf, 1.0)
                fv <- ores$fval
                sig <- sqrt((rr$wrss() + pp1$sqrL(1))/n)
                c(sign(fw - est) * sqrt(fv - base),
                  Cv_to_Sv(ores$par, vapply(fitted@cnms,length, 1), s=sig),
                  ## ores$par * sig, sig,
                  mkpar(p, j, fw, pp1$beta(1)))
            }
            nres[1, ] <- pres[2, ] <- fe.zeta(est + delta * std)
            poff <- nvp + 1L + j
            bres <-
                as.data.frame(unique(rbind2(fillmat(pres,-Inf, Inf, fe.zeta, poff),
                                            fillmat(nres,-Inf, Inf, fe.zeta, poff))))
            thisnm <- names(fe.orig)[j]
            bres$.par <- thisnm
            ans[[thisnm]] <- bres[order(bres[, poff]), ]
            form[[3]] <- as.name(thisnm)
            bakspl[[thisnm]] <-
                tryCatch(backSpline(forspl[[thisnm]] <- interpSpline(form, bres)),
                         error=function(e)e)
            if (inherits(bakspl[[thisnm]],"error")) warning("non-monotonic profile")
        } ## for(j in 1..p)
    } ## if isLMM

    ans <- do.call(rbind, ans)
    row.names(ans) <- NULL
    ans$.par <- factor(ans$.par)
    attr(ans, "forward") <- forspl
    attr(ans, "backward") <- bakspl
    class(ans) <- c("thpr", "data.frame")
    attr(ans, "lower") <- lower[seqnvp]
    attr(ans, "upper") <- upper[seqnvp]
    ans
} ## profile.merMod

get.which <- function(which,nvp,nptot,parnames,verbose) {
    wi.vp <- seq_len(nvp)
    if(is.character(which)) {
	wi <- integer(); wh <- which
	if(any(j <- wh == "theta_")) {
	    wi <- wi.vp; wh <- wh[!j] }
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
        which <- wi
    } # else stopifnot( .. numeric ..)
    which
}

## This is a hack.  The preferred approach is to write a
## subset method for the ddenseModelMatrix and dsparseModelMatrix
## classes
.modelMatrixDrop <- function(mm, w) {
    if (isS4(mm)) {
        ll <- list(Class = class(mm),
                   assign = attr(mm,"assign")[-w],
                   contrasts = NULL)
        ## FIXME: where did the contrasts information go??
        ## mm@contrasts)
        X <- mm[, -w, drop = FALSE]
        ll <- c(ll, lapply(structure(slotNames(X), .Names=slotNames(X)),
                           function(nm) slot(X, nm)))
        return(do.call("new", ll))
    }
    ans <- mm[, -w, drop=FALSE]
    attr(ans, "assign") <- attr(mm, "assign")[-w]
    ans
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
    basedev <- deviance(fm)
    vlist <- sapply(fm@cnms,length)
    sig <- sigma(fm)  ## only if useSc=TRUE?
    stdErr <- unname(coef(summary(fm))[,2])
    pp <- fm@pp$copy()
    ## opt <- c(pp$theta*sig, sig)
    if (useSc) {
        opt <- Cv_to_Sv(pp$theta, n=vlist, s=sig)
        names(opt) <- if (signames) {
            c(sprintf(".sig%02d", seq(length(opt)-1)), ".sigma")
        } else {
            c(tnames(fm,old=FALSE,prefix=c("sd","cor")),"sigma")
        }
    } else {
        opt <- Cv_to_Sv(pp$theta, n=vlist)
        names(opt) <- if (signames) {
            sprintf(".sig%02d", seq_along(opt))
        } else {
            tnames(fm,old=FALSE,prefix=c("sd","cor"))
        }
    }
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
        panel.lines(myspl$x,myspl$y)
        panel.points(myspl$x,myspl$y,pch="+")
        warning(sprintf("bad profile for variable %d: using linear interpolation",panel.number()))
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
              levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 1)),
              conf = c(50, 80, 90, 95, 99)/100,
              absVal = FALSE,
              which = 1:nptot, ...)
{
    nptot <- length(nms <- levels(x[[".par"]]))
    ## FIXME: is this sufficiently reliable?
    ## (include "sigma" in 'theta' parameters)
    nvp <- length(grep("^(\\.sig[0-9]+|.sigma|sd_|cor_)",nms))
    which <- get.which(which,nvp,nptot,nms,verbose=FALSE)
    levels <- sort(levels[is.finite(levels) & levels > 0])
    spl <- attr(x, "forward")[which]
    bspl <- attr(x, "backward")[which]
    ## for parameters for which spline couldn't be computed,
    ## replace the 'spl' element with the raw profile data
    badSpl <- sapply(spl,is.null)
    spl[badSpl] <- lapply(which(badSpl),
                          function(i) {
                              n <- names(badSpl)[i]
                              r <- x[x[[".par"]]==n,]
                              data.frame(y=r[[".zeta"]],
                                         x=r[[n]])
                          })
    bspl[badSpl] <- lapply(spl[badSpl],function(d) { data.frame(x=d$y,y=d$x) })
    zeta <- c(-rev(levels), 0, levels)
    ## use linear approximation if backspline doesn't work
    tmpf <- function(bspl,zeta) {
        if (is(bspl,"spline")) return(predy(bspl,zeta))
        return(approx(bspl$x,bspl$y,xout=zeta)$y)
    }
    fr <- data.frame(zeta = rep.int(zeta, length(spl)),
                     pval = unlist(lapply(bspl,tmpf,zeta)),
                     pnm = gl(length(spl), length(zeta), labels = names(spl)))
    if (length(ind <- which(is.na(fr$pval)))) {
        fr[ind, "zeta"] <- 0
        for (i in ind)
### FIXME: Should check which bound has been violated, although it
### will almost always be the minimum.
            if (is(curspl <-  spl[[fr[i, "pnm"] ]],"spline")) {
                fr[i, "pval"] <- min(curspl$knots)
            }
    }
    ylab <- expression(zeta)
    if (absVal) {
        fr$zeta <- abs(fr$zeta)
        ylab <- expression("|" * zeta * "|")
    }
    ll <- c(list(...),
            list(x = zeta ~ pval | pnm, data=fr,
                 scales = list(x = list(relation = 'free')),
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

##' @importFrom stats confint
##' @S3method confint thpr
confint.thpr <- function(object, parm, level = 0.95, zeta, ...)
{
    bak <- attr(object, "backward")
    ## fallback strategy for old profiles that don't have a lower/upper
    ##  attribute saved ...
    if (is.null(lower <- attr(object,"lower")))
        lower <- rep(NA,length(parm))
    if (is.null(upper <- attr(object,"upper")))
        upper <- rep(NA,length(parm))
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
        nm <- parm[i]
        ## would like to build this machinery into predy, but
        ## predy is used in many places and it's much harder to
        ## tell in general whether an NA indicates a lower or an
        ## upper bound ...
        if (!is(b <- bak[[nm]],"try-error")) {
            p <- predy(b, zeta)
            if (is.na(p[1])) p[1] <- lower[i]
            if (is.na(p[2])) p[2] <- upper[i]
            ci[i,] <- p
        }
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

##' @details Depending on the method specified, this function will
##' compute confidence intervals by ("profile") computing a likelihood profile
##' and finding the appropriate cutoffs based on the likelihood ratio test;
##' ("Wald") approximate the confidence intervals (of fixed-effect parameters
##' only) based on the estimated local curvature of the likelihood surface;
##' ("boot") perform parametric bootstrapping
##' with confidence intervals computed from the bootstrap distribution
##' according to \code{boot.type} (see \code{\link{boot.ci}})
##' @importFrom stats confint
##' @S3method confint merMod
##' @method confint merMod
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm1W <- confint(fm1,method="Wald")
##' \dontrun{
##' ## ~20 seconds, MacBook Pro laptop
##' system.time(fm1P <- confint(fm1,method="profile",oldNames=FALSE)) ## default
##' ## ~ 40 seconds
##' system.time(fm1B <- confint(fm1,method="boot",
##'                     .progress="txt", PBargs=list(style=3)))
##' }
##' load(system.file("testdata","confint_ex.rda",package="lme4"))
##' fm1P
##' fm1B
confint.merMod <- function(object, parm, level = 0.95,
			   method=c("profile","Wald","boot"),
			   zeta, nsim=500, boot.type=c("perc","basic","norm"),
                           quiet=FALSE, oldNames=TRUE, ...)
{
    method <- match.arg(method)
    boot.type <- match.arg(boot.type)
    if (!missing(parm) && !is.numeric(parm) && method %in% c("profile","boot"))
        stop("for method='",method,"', 'parm' must be specified as an integer")
    switch(method,
	   "profile" =
       {
           if (!quiet) message("Computing profile confidence intervals ...")
           utils::flush.console()
	   pp <- if(missing(parm)) {
               profile(object, signames=oldNames, ...)
           } else {
               profile(object, which=parm, signames=oldNames, ...)
           }
	   confint(pp,level=level,zeta=zeta)
       },
	   "Wald" =
       {
        ## copied with small changes from confint.default
        cf <- fixef(object)   ## coef() -> fixef()
        pnames <- names(cf)
        if (missing(parm))
          parm <- pnames
        else if (is.numeric(parm))
          parm <- pnames[parm]
        ## n.b. can't use sqrt(...)[parm] (diag() loses names)
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        pct <- format.perc(a, 3)
        fac <- qnorm(a)
        ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                                    pct))
        sdiag <- function(x) if (length(x)==1) c(x) else diag(x)
        ses <- sqrt(sdiag(vcov(object)[parm,parm]))
        ci[] <- cf[parm] + ses %o% fac
        ci
        ## only gives confidence intervals on fixed effects ...
    },
	   "boot" =
       {
	   if (!quiet) message("Computing bootstrap confidence intervals ...")
           utils::flush.console()
           bootFun <- function(x) {
               th <- getME(x,"theta")
               nvec <- sapply(getME(x,"cnms"),length)
               scaleTh <- (isLMM(x) || isNLMM(x))
               useSc <- as.logical(x@devcomp$dims["useSc"])
               ## FIXME: still ugly.  Best cleanup via Cv_to_Sv ...
               tn <- tnames(x,old=FALSE,prefix=c("sd","cor"))
               if (scaleTh) {  ## scale variances by sigma and include it
                   ss <- setNames(Cv_to_Sv(th,n=nvec,s=sigma(x)),
                                  c(tn,"sigma"))
               } else if (useSc) { ## don't scale variances but do include sigma
                   ss <- setNames(c(Cv_to_Sv(th,n=nvec),sigma(x)),
                                  c(tn,"sigma"))
               } else {  ## no scaling, no sigma
                   ss <- setNames(Cv_to_Sv(th,n=nvec),
                                  tn)
               }
               c(ss, fixef(x))
           }
	   bb <- bootMer(object, bootFun, nsim=nsim,...)
           bci <- lapply(seq_along(bb$t0),
                         boot.out=bb,
                         boot::boot.ci,type=boot.type,conf=level)
           cpos <- grep(boot.type,names(bci[[1]]))
           ## get _last_ two columns
           ccol <- ncol(bci[[1]][[cpos]])+(-1:0)
           citab <- t(sapply(bci,function(x) x[[cpos]][ccol]))
           a <- (1 - level)/2
           a <- c(a, 1 - a)
           pct <- format.perc(a, 3)
           dimnames(citab) <- list(names(bb[["t0"]]),pct)
           pnames <- rownames(citab)
           if (missing(parm))
               parm <- pnames
           else if (is.numeric(parm))
               parm <- pnames[parm]
           citab[parm,]
       },
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
                    conf = c(50, 80, 90, 95, 99)/100,
                    which=1:nptot, ...)
{
    singfit <- FALSE
    for (i in grep("^(\\.sig[0-9]+|sd_)",names(x))) {
        singfit <- singfit ||  (any(x[,".zeta"]==0 & x[,i]==0))
    }
    if (singfit) warning("splom is unreliable for singular fits")

    nptot <- length(nms <- names(attr(x,"forward")))
    nvp <- length(grep("^(\\.sig[0-9]+|.sigma|sd_|cor_)",nms))
    which <- get.which(which,nvp,nptot,nms,verbose=FALSE)

    if (length(which)==1)
        stop("can't draw a scatterplot matrix for a single variable")

    mlev <- max(levels)
    spl <- attr(x, "forward")[which]
    frange <- sapply(spl, function(x) range(x$knots))
    bsp <- attr(x, "backward")[which]
    x <- x[x[[".par"]] %in% nms[which],c(".zeta",nms[which],".par")]
    brange <- sapply(bsp, function(x) range(x$knots))
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
    ## panel function for lower triangle
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
    ## panel function for upper triangle
    up <- function(x, y, groups, subscripts, i, j, ...) {
        ## panels are transposed so reverse i and j
        jj <- i
        ii <- j
        tr <- traces[[jj]][[ii]]
        ll <- tr$ll
        pts <- ll$pts
        limits <- current.panel.limits()
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

    splom(~ pfr, lower.panel = lp, upper.panel = up, diag.panel = dp, ...)
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


