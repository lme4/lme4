##' Methods for profile() of [ng]lmer fitted models
##' 
##' Methods for function \code{\link{profile}} (package \pkg{stats}), here for
##' profiling (fitted) mixed effect models.
##' 
##' 
##' @name profile-methods
##' @aliases profile-methods profile.merMod
##' @docType methods
##' @param fitted a fitted model, e.g., the result of \code{\link{lmer}(..)}.
##' @param alphamax used when \code{delta} is unspecified, as probability ... to
##' compute \code{delta} ...
##' @param maxpts ...
##' @param delta ...
##' @param tr ...
##' @param \dots potential further arguments for \code{profile} methods.
##' @section Methods: FIXME: These (signatures) will change soon --- document
##' \bold{after} change!
##' \describe{
##'     \item{signature(fitted = \"merMod\")}{ ...  } }
##' @seealso For (more expensive) alternative confidence intervals:
##' \code{\link{bootMer}}.
##' @keywords methods
##' @examples
##' fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
##' ## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
##' system.time( tpr <- profile(fm01ML) )
##' (confint(tpr) -> CIpr)
##' print(xyplot(tpr))
##' @importFrom splines backSpline interpSpline periodicSpline
##' @importFrom stats profile
##' @method profile merMod
##' @export
profile.merMod <- function(fitted, alphamax = 0.01, maxpts = 100, delta = cutoff/8,
                           tr = 0, ...) {
    dd <- devfun2(fitted)
    
    base <- attr(dd, "basedev")
    thopt <- attr(dd, "thopt")
    stderr <- attr(dd, "stderr")
    pp <- environment(dd)$pp
    X.orig <- pp$X
    n <- environment(dd)$n
    p <- length(pp$beta0)
    
    ans <- lapply(opt <- attr(dd, "optimum"), function(el) NULL)
    bakspl <- forspl <- ans
    
    nptot <- length(opt)
    nvp <- nptot - p    # number of variance-covariance pars
    fe.orig <- opt[-seq_len(nvp)]
    res <- c(.zeta = 0, opt)
    res <- matrix(res, nr = maxpts, nc = length(res),
                  dimnames = list(NULL, names(res)), byrow = TRUE)
    cutoff <- sqrt(qchisq(1 - alphamax, nptot))
    
    ## helper functions
    
    ## nextpar calculates the next value of the parameter being
    ## profiled based on the desired step in the profile zeta
    ## (absstep) and the values of zeta and column cc for rows
    ## r-1 and r.  The parameter may not be below lower
    nextpar <- function(mat, cc, r, absstep, lower = -Inf) {
        rows <- r - (1:0)         # previous two row numbers
        pvals <- mat[rows, cc]
        zeta <- mat[rows, ".zeta"]
        if (!(denom <- diff(zeta)))
            stop("Last two rows have identical .zeta values")
        num <- diff(pvals)
        max(lower, pvals[2] + sign(num) * absstep * num / denom)
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
    fillmat <- function(mat, lowcut, zetafun, cc) {
        nr <- nrow(mat)
        i <- 2L
        while (i < nr && abs(mat[i, ".zeta"]) <= cutoff &&
               mat[i, cc] > lowcut) {
            mat[i + 1L, ] <- zetafun(nextpar(mat, cc, i, delta, lowcut))
            i <- i + 1L
        }
        mat
    }
    
    lower <- c(fitted@lower, 0, rep.int(-Inf, p))
    seqnvp <- seq_len(nvp)
    lowvp <- lower[seqnvp]
    form <- .zeta ~ foo           # pattern for interpSpline formula
    
    for (w in seqnvp) {
        wp1 <- w + 1L
        start <- opt[seqnvp][-w]
        pw <- opt[w]
        lowcut <- lower[w]
        zeta <- function(xx) {
            ores <- bobyqa(start,
                           function(x) dd(mkpar(nvp, w, xx, x)),
                           lower = lowvp[-w])
            zz <- sign(xx - pw) * sqrt(ores$fval - base)
            c(zz, mkpar(nvp, w, xx, ores$par), pp$beta(1))
        }
        
### FIXME: The starting values for the conditional optimization should
### be determined from recent starting values, not always the global
### optimum values.
        
### Can do this by taking the change in the other parameter values at
### the two last points and extrapolating.
        
        ## intermediate storage for pos. and neg. increments
        pres <- nres <- res
        ## assign one row, determined by inc. sign, from a small shift
        nres[1, ] <- pres[2, ] <- zeta(pw * 1.01)
        ## fill in the rest of the arrays and collapse them
        bres <-
            as.data.frame(unique(rbind2(fillmat(pres,lowcut, zeta, wp1),
                                        fillmat(nres,lowcut, zeta, wp1))))
        bres$.par <- names(opt)[w]
        ans[[w]] <- bres[order(bres[, wp1]), ]
        form[[3]] <- as.name(names(opt)[w])
        
        bakspl[[w]] <- backSpline(forspl[[w]] <- interpSpline(form, bres))
    }

    offset.orig <- fitted@resp$offset
    for (j in seq_len(p)) {
        pres <-            # intermediate results for pos. incr.
            nres <- res    # and negative increments
        est <- opt[nvp + j]
        std <- stderr[j]
        Xw <-X.orig[, j, drop=TRUE]
        Xdrop <- .modelMatrixDrop(X.orig, j)
        pp1 <- do.call(new, list(Class = class(pp),
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
        fe.zeta <- function(fw) {
            rr$setOffset(Xw * fw + offset.orig)
            rho <- as.environment(list(pp=pp1, resp=rr))
            parent.env(rho) <- parent.frame()
            ores <- bobyqa(thopt, mkdevfun(rho, 0L), lower = fitted@lower)
            fv <- ores$fval
            sig <- sqrt((rr$wrss() + pp1$sqrL(1))/n)
            c(sign(fw - est) * sqrt(fv - base),
              ores$par * sig, sig, mkpar(p, j, fw, pp1$beta(1)))
        }
        nres[1, ] <- pres[2, ] <- fe.zeta(est + delta * std)
        poff <- nvp + 1L + j
        bres <-
            as.data.frame(unique(rbind2(fillmat(pres,-Inf, fe.zeta, poff),
                                        fillmat(nres,-Inf, fe.zeta, poff))))
        thisnm <- names(fe.orig)[j]
        bres$.par <- thisnm
        ans[[thisnm]] <- bres[order(bres[, poff]), ]
        form[[3]] <- as.name(thisnm)
        bakspl[[thisnm]] <-
            backSpline(forspl[[thisnm]] <- interpSpline(form, bres))
    }
    
    ans <- do.call(rbind, ans)
    ans$.par <- factor(ans$.par)
    attr(ans, "forward") <- forspl
    attr(ans, "backward") <- bakspl
    row.names(ans) <- NULL
    class(ans) <- c("thpr", "data.frame")
    ans
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
devfun2 <- function(fm)
{
    stopifnot(is(fm, "lmerMod"))
    fm <- refitML(fm)
    basedev <- deviance(fm)
    sig <- sigma(fm)
    stdErr <- unname(coef(summary(fm))[,2])
    pp <- fm@pp$copy()
    opt <- c(sig * pp$theta, sig)
    names(opt) <- c(sprintf(".sig%02d", seq_along(pp$theta)), ".sigma")
    opt <- c(opt, fixef(fm))
    resp <- fm@resp$copy()
    np <- length(pp$theta) + 1L
    n <- nrow(pp$V)                   # use V, not X so it works with nlmer
    ans <- function(pars)
    {
        stopifnot(is.numeric(pars), length(pars) == np)
        ## Assumption:  all parameters, including the residual SD on SD-scale
        sigma <- pars[np]
        .Call(lme4Eigen:::lmer_Deviance, pp$ptr(), resp$ptr(), pars[-np]/sigma)
        sigsq <- sigma^2
        pp$ldL2() + (resp$wrss() + pp$sqrL(1))/sigsq + n * log(2 * pi * sigsq)
    }
    attr(ans, "optimum") <- opt         # w/ names()
    attr(ans, "basedev") <- basedev
    attr(ans, "thopt") <- pp$theta
    attr(ans, "stderr") <- stdErr
    class(ans) <- "devfun"
    ans
}

## extract only the y component from a prediction
predy <- function(sp, vv) predict(sp, vv)$y

stripExpr <- function(ll, nms) {
    stopifnot(inherits(ll, "list"), is.character(nms))
    sigNm <- which(nms == ".sigma")
    sigNms <- grep("^.sig[0-9]+", nms)
    sigsub <- as.integer(substring(nms[sigNms], 5))
    fLevs <- as.expression(nms)
    fLevs[sigNm] <- expression(sigma)
    fLevs[sigNms] <- parse(text=paste("sigma[", sigsub, "]"))
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

## A lattice-based plot method for profile objects
##' @importFrom lattice xyplot
##' @S3method xyplot thpr
xyplot.thpr <-
    function (x, data = NULL,
              levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 1)),
              conf = c(50, 80, 90, 95, 99)/100,
              absVal = FALSE, ...)
{
    levels <- sort(levels[is.finite(levels) & levels > 0])
    spl <- attr(x, "forward")
    bspl <- attr(x, "backward")
    zeta <- c(-rev(levels), 0, levels)
    fr <- data.frame(zeta = rep.int(zeta, length(spl)),
                     pval = unlist(lapply(bspl, predy, zeta)),
                     pnm = gl(length(spl), length(zeta), labels = names(spl)))
    if (length(ind <- which(is.na(fr$pval)))) {
        fr[ind, "zeta"] <- 0
        for (i in ind)
### FIXME: Should check which bound has been violated, although it
### will almost always be the minimum.
            fr[i, "pval"] <- min(spl[[fr[i, "pnm"] ]]$knots)
    }
    ylab <- expression(zeta)
    if (absVal) {
        fr$zeta <- abs(fr$zeta)
        ylab <- expression("|" * zeta * "|")
    }
    ll <- c(list(...),
            list(x = zeta ~ pval | pnm, data=fr,
                 scales = list(x = list(relation = 'free')),
                 ylab = ylab, xlab = NULL,
                 panel = function(x, y, ...)
             {
                 panel.grid(h = -1, v = -1)
                 lsegments(x, y, x, 0, ...)
                 lims <- current.panel.limits()$xlim
                 myspl <- spl[[panel.number()]]
                 krange <- range(myspl$knots)
                 pr <- predict(myspl,
                               seq(max(lims[1], krange[1]),
                                   min(lims[2], krange[2]), len = 101))
                 if (absVal) {
                     pr$y <- abs(pr$y)
                     y[y == 0] <- NA
                     lsegments(x, y, rev(x), y)
                 } else {
                     panel.abline(h = 0, ...)
                 }
                 panel.lines(pr$x, pr$y)
             }))
    do.call(xyplot, stripExpr(ll, names(spl)))
}

##' @importFrom stats confint
##' @S3method confint thpr
confint.thpr <- function(object, parm, level = 0.95, zeta, ...)
{
    bak <- attr(object, "backward")
    bnms <- names(bak)
    if (missing(parm)) parm <- bnms
    else if (is.numeric(parm)) parm <- bnms[parm]
    parm <- intersect(as.character(parm), bnms)
    cn <- NULL
    if (missing(zeta)) {
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        zeta <- qnorm(a)
        cn <- stats:::format.perc(a, 3)
    }
    ci <- t(sapply(parm, function(nm) predy(bak[[nm]], zeta)))
    colnames(ci) <- cn
    ci
}

## Convert x-cosine and y-cosine to average and difference.

## Convert the x-cosine and the y-cosine to an average and difference
## ensuring that the difference is positive by flipping signs if
## necessary
## @param xc x-cosine
## @param yc y-cosine
ad <- function(xc, yc)
{
    a <- (xc + yc)/2
    d <- (xc - yc)
    cbind(ifelse(d > 0, a, -a), abs(d))
}

## convert d versus a (as an xyVector) and level to a matrix of taui and tauj
## @param xy an xyVector
## @param lev the level of the contour
tauij <- function(xy, lev) lev * cos(xy$x + outer(xy$y/2, c(-1, 1)))

## @title safe arc-cosine
## @param x numeric vector argument
## @return acos(x) being careful of boundary conditions
sacos <- function(x) acos(pmax.int(-0.999, pmin.int(0.999, x)))

## Generate a contour
##
## @title Generate a contour
## @param sij the arc-cosines of i on j
## @param sji the arc-cosines of j on i
## @param levels numeric vector of levels at which to interpolate
## @param nseg number of segments in the interpolated contour
## @return a list with components
## \item{tki}{the tau-scale predictions of i on j at the contour levels}
## \item{tkj}{the tau-scale predictions of j on i at the contour levels}
## \item{pts}{an array of dimension (length(levels), nseg, 2) containing the points on the contours}
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


## Profile pairs plot
## Contours are for the marginal two-dimensional regions (i.e. using
##  df = 2)
##' @importFrom grid gpar viewport
##' @importFrom lattice splom
##' @S3method splom thpr
splom.thpr <-
    function (x, data, ## unused - only for compatibility with generic
              levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 2)),
              conf = c(50, 80, 90, 95, 99)/100, ...)
{
    mlev <- max(levels)
    spl <- attr(x, "forward")
    frange <- sapply(spl, function(x) range(x$knots))
    bsp <- attr(x, "backward")
    brange <- sapply(bsp, function(x) range(x$knots))
    pfr <- do.call(cbind, lapply(bsp, predy, c(-mlev, mlev)))
    pfr[1, ] <- pmax.int(pfr[1, ], frange[1, ], na.rm = TRUE)
    pfr[2, ] <- pmin.int(pfr[2, ], frange[2, ], na.rm = TRUE)
    nms <- names(spl)
    ## Create data frame fr of par. vals in zeta coordinates
    fr <- x[, -1]
    for (nm in nms) fr[[nm]] <- predy(spl[[nm]], fr[[nm]])
    fr1 <- fr[1, nms]
    ## create a list of lists with the names of the parameters
    traces <- lapply(fr1, function(el) lapply(fr1, function(el1) list()))
    for (j in seq_along(nms)[-1]) {
        for (i in seq_len(j - 1)) {
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
                           fontface = lattice:::chooseFace(varname.fontface,
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
                           tick = TRUE,
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
            side <- ifelse(j == 1, "right", "bottom")
            which.half <- ifelse(j == 1, "lower", "upper")
            at <- pretty(lims)
            panel.axis(side = side, at = at, labels = format(at, trim = TRUE),
                       tick = TRUE, half = TRUE, which.half = which.half,
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

## Transform an lmer profile to the scale of the logarithm of the
## standard deviation of the random effects.
## @title Transform an lmer profile to the logarithm scale
## @param x an object that inherits from class "thpr"
## @param base the base of the logarithm.  Defaults to natural
##        logarithms
##
## @return an lmer profile like x with all the .sigNN parameters
##      replaced by .lsigNN.  The forward and backward splines for
##      these parameters are recalculated.
##' @S3method log thpr
log.thpr <- function (x, base = exp(1)) {
    cn <- colnames(x)
    sigs <- grep("^\\.sig", cn)
    if (length(sigs)) {
        colnames(x) <- sub("^\\.sig", ".lsig", cn)
        levels(x$.par) <- sub("^\\.sig", ".lsig", levels(x$.par))
        names(attr(x, "backward")) <-
            names(attr(x, "forward")) <-
                sub("^\\.sig", ".lsig", names(attr(x, "forward")))
        for (nm in colnames(x)[sigs]) {
            x[[nm]] <- log(x[[nm]], base = base)
            fr <- subset(x, .par == nm & is.finite(x[[nm]]))
            ## FIXME: avoid subset for global-variable false positive
            ## fr <- x[x$.par == nm & is.finite(x[[nm]]),]
            form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
            attr(x, "forward")[[nm]] <- interpSpline(form, fr)
            attr(x, "backward")[[nm]] <- backSpline(attr(x, "forward")[[nm]])
        }
        ## eliminate rows the produced non-finite logs
        x <- x[apply(is.finite(as.matrix(x[, sigs])), 1, all),]
    }
    x
}

## Transform a profile from the standard deviation parameters to the variance
##
## @title Transform to variance component scale
## @param x a profile object from a mixed-effects model
## @return a modified profile object
varpr <- function (x) {
    cn <- colnames(x)
    sigs <- grep("^\\.sig", cn)
    if (length(sigs)) {
        colnames(x) <- sub("^\\.sig", ".sigsq", cn)
        levels(x$.par) <- sub("^\\.sig", ".sigsq", levels(x$.par))
        names(attr(x, "backward")) <-
            names(attr(x, "forward")) <-
                sub("^\\.sig", ".sigsq", names(attr(x, "forward")))
        for (nm in colnames(x)[sigs]) {
            x[[nm]] <- x[[nm]]^2
            fr <- subset(x, .par == nm & is.finite(x[[nm]]))
            form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
            attr(x, "forward")[[nm]] <- interpSpline(form, fr)
            attr(x, "backward")[[nm]] <- backSpline(attr(x, "forward")[[nm]])
        }
        ## eliminate rows the produced non-finite logs
        x <- x[apply(is.finite(as.matrix(x[, sigs])), 1, all),]
    }
    x
}

## Create an approximating density from a profile object
##
## @title Approximate densities from profiles
## @param pr a profile object
## @param npts number of points at which to evaluate the density
## @param upper upper bound on cumulative for a cutoff
## @return a data frame
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

## Densityplot method for a mixed-effects model profile
##
## @title densityplot from a mixed-effects profile
## @param x a mixed-effects profile
## @param data not used - for compatibility with generic
## @param ... optional arguments to densityplot
## @return a density plot
##' @importFrom lattice densityplot
##' @S3method densityplot thpr
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
##'
##' @title Transform to the variance scale
##' @param pr a mixed-effects model profile
##' @return a transformed mixed-effects model profile
##' @export
varianceProf <- function(pr) {
    stopifnot(inherits(pr, "thpr"))
    spl <- attr(pr, "forward")
    onms <- names(spl)                  # names of original variables
    vc <- onms[grep("^.sig", onms)]     # variance components
    ans <- subset(pr, .par %in% vc, select=c(".zeta", vc, ".par"))
    ans$.par <- factor(ans$.par)        # drop unused levels
    if (".lsig" %in% vc) ans$.lsig <- exp(ans$.lsig)
    attr(ans, "forward") <- attr(ans, "backward") <- list()
    for (nm in vc) {
        ans[[nm]] <- ans[[nm]]^2
        fr <- subset(ans, .par == nm & is.finite(ans[[nm]]))
        form <- eval(substitute(.zeta ~ nm, list(nm = as.name(nm))))
        attr(ans, "forward")[[nm]] <- interpSpline(form, fr)
        attr(ans, "backward")[[nm]] <- backSpline(attr(ans, "forward")[[nm]])
    }
    ans
}
