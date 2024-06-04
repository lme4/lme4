### Plots for the ranef.mer class ----------------------------------------

##' @importFrom lattice dotplot
##' @S3method  dotplot ranef.mer
dotplot.ranef.mer <- function(x, data, main = TRUE, transf=I, level = 0.95, ...)
{
    rng <- qnorm((1+level)/2)
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        if (is.null(se)) return(list())
        x <- as.numeric(x)
        hw <- rng * as.numeric(se[subscripts])
        list(xlim = range(transf(x - hw), transf(x + hw), finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16,
                         horizontal = TRUE, col = dot.symbol$col, 
                         lty.h = dot.line$lty, lty.v = dot.line$lty, 
                         lwd.h = dot.line$lwd, lwd.v = dot.line$lwd, 
                         col.line.h = dot.line$col, col.line.v = dot.line$col, 
                         levels.fos = unique(y),
                         groups = NULL, ...)
    {
        x <- as.numeric(x)
        y <- as.numeric(y)
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        sup.symbol <- trellis.par.get("superpose.symbol")
        panel.abline(h = levels.fos, col = col.line.h, lty = lty.h, lwd = lwd.h)
        panel.abline(v = 0, col = col.line.v, lty = lty.v, lwd = lwd.v)
        if (!is.null(se)) {
            se <- as.numeric(se[subscripts])
            panel.segments( transf(x - rng * se), y,
                            transf(x + rng * se), y, col = 'black')
        }
        panel.xyplot(transf(x), y, pch = pch, col = col, ...)
    }
    f <- function(nx, ...) {
        ss <- asDf0(x,nx)
        mtit <- if(main) nx 
        dotplot(.nn ~ values | ind, ss, se = ss$se,
                prepanel = prepanel.ci, panel = panel.ci,
                xlab = NULL, main = mtit, ...)
    }
    setNames(lapply(names(x), f, ...), names(x))
}

##' @importFrom graphics plot
##' @S3method plot ranef.mer
plot.ranef.mer <- function(x, y, ...)
{
    lapply(x, function(x) {
        cn <- lapply(colnames(x), as.name)
        switch(min(ncol(x), 3),
               qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
               xyplot(eval(substitute(y ~ x,
                                      list(y = cn[[1]],
                                           x = cn[[2]]))), x, ...),
               splom(~ x, ...))
    })
}

##' @importFrom lattice qqmath
##' @S3method qqmath ranef.mer
qqmath.ranef.mer <- function(x, data, main = TRUE, level = 0.95, ...)
{
    rng <- qnorm((1+level)/2)
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        x <- as.numeric(x)
        se <- as.numeric(se[subscripts])
        hw <- rng * se
        list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
        panel.grid(h = -1,v = -1)
        panel.abline(v = 0)
        x <- as.numeric(x)
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        panel.segments(x - rng * se, y, x + rng * se, y, col = 'black')
        panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(nx) {
        xt <- x[[nx]]
        mtit <- if(main) nx # else NULL
        if (!is.null(pv <- attr(xt, "postVar")))
        {
            d <- dim(pv)
            se <- vapply(seq_len(d[1]), function(i) sqrt(pv[i, i, ]), numeric(d[3]))
            nr <- nrow(xt)
            nc <- ncol(xt)
            ord <- unlist(lapply(xt, order)) + rep((0:(nc - 1)) * nr, each = nr)
            rr <- 1:nr
            ind <- gl(nc, nr, labels = names(xt))
            xyplot(rep(qnorm((rr - 0.5)/nr), nc) ~ unlist(xt)[ord] | ind[ord],
                   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                   scales = list(x = list(relation = "free")),
                   ylab = "Standard normal quantiles",
                   xlab = NULL, main = mtit, ...)
        } else {
            qqmath(~values|ind,
                   data = stack(xt),
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, main = mtit, ...)
        }
    }
    sapply(names(x), f, simplify = FALSE)
}

##' @importFrom graphics plot
##' @S3method plot coef.mer
plot.coef.mer <- function(x, y, ...)
{
    ## remove non-varying columns from frames
    reduced <- lapply(x, function(el)
                      el[, !vapply(el, function(cc) all(cc == cc[1L]), NA)])
    plot.ranef.mer(reduced, ...)
}

##' @importFrom lattice dotplot
##' @S3method dotplot coef.mer
dotplot.coef.mer <- function(x, data, ...) {
    mc <- match.call()
    mc[[1]] <- as.name("dotplot.ranef.mer")
    eval(mc)
}
