## List of linear models according to a grouping factor

## Extract the model formula
modelFormula <- function(form)
{
    if (!inherits(form, "formula") || length(form) != 3)
        stop("formula must be a two-sided formula object")
    rhs <- form[[3]]
    if (!inherits(rhs, "call") || rhs[[1]] != as.symbol('|'))
        stop("rhs of formula must be a conditioning expression")
    form[[3]] <- rhs[[2]]
    list(model = dropOffset(form), groups = rhs[[3]])
}

dropOffset <- function(form) {
    ## atomic
    if (is.symbol(form) || is.numeric(form)) return(form)
    ## binary
    if (identical(form[[1]],quote(offset))) {
        NULL
    } else {
        ## unary operator
        if (length(form)==2) {
            form[[2]] <- dropOffset(form[[2]])
            return(form)
        }
        nb2 <- dropOffset(form[[2]])
        nb3 <- dropOffset(form[[3]])
        if (is.null(nb2))
            nb3
        else if (is.null(nb3))
            nb2
        else {
            form[[2]] <- nb2
            form[[3]] <- nb3
            return(form)
        }
    }
}

## dropOffset(y~x)
## dropOffset(y~x+offset(stuff))
## dropOffset(y~-x+offset(stuff))
## dropOffset(~-x+offset(stuff))

##' @title List of lm Objects with a Common Model
##' @param formula a linear formula object of the form
##'     \code{y ~ x1+...+xn | g}. In the formula object, \code{y} represents
##'     the response, \code{x1,...,xn} the covariates, and \code{g} the
##'     grouping factor specifying the partitioning of the data according to
##'     which different \code{lm} fits should be performed.
##' @inheritParams lmer
##' @param family an optional family specification for a generalized
##'     linear model.
##' @param pool logical scalar, should the variance estimate pool the
##'     residual sums of squares
##' @param ... additional, optional arguments to be passed to the
##'     model function or family evaluation.
##' @export
lmList <- function(formula, data, family, subset, weights,
                   na.action, offset, pool = TRUE, ...) {
    stopifnot(inherits(formula, "formula"))

    ## model.frame(groupedData) was problematic ... but not as we
    ## are currently using it.

    mCall <- mf <- match.call()

    ## in contrast to the usual R model-fitting idiom, we do **not**
    ## want to evaluate the model frame here; it will mess up any derived
    ## variables we have when we go to fit the sub-models.  We were previously
    ## using model.frame() on the entire data set, but that does not
    ## exclude unused columns ... and hence screws us up when there are
    ## NA values in unused columns.  All we need the model frame for
    ## is evaluating the groups.

    ## keep weights and offsets in case we have NAs there??
    m <- match(c("formula", "data", "subset", "na.action",
                 "weights", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    ## substitute `+' for `|' in the formula
    mf$formula <- subbars(formula)
    mf$drop.unused.levels <- TRUE
    ## pass NAs for now -- want *all* groups, weights, offsets recovered
    mf$na.action <- na.pass
    mf[[1L]] <- quote(stats::model.frame)
    frm <- eval.parent(mf) ## <- including "..."
    data[["(weights)"]] <- model.weights(frm)
    data[["(offset)"]] <- model.offset(frm)
    mform <- modelFormula(formula)
    isGLM <- !missing(family) ## TODO in future, consider isNLM / isNLS
    errorH <- function(e) NULL # => NULL iff an error happened
    ## FIXME: catch errors and pass them on as warnings?
    ## (simply passing them along silently gives confusing output)
    groups <- eval(mform$groups, frm)
    if (!is.factor(groups)) groups <- factor(groups)
    fit <- if (isGLM) glm else lm
    mf2 <- if (missing(family)) NULL else list(family=family)
    fitfun <- function(data,formula) {
        tryCatch({
            do.call(fit,c(list(formula, data,
                               weights = model.weights(data),
                               offset = model.offset(data), ...),
                          mf2))
        }, error=errorH)
    }
    ## split *original data*, not model frame, on groups
    frm.split <- split(data, groups)
    ## NB:  levels() is only  OK if grouping variable is a factor
    nms <- names(frm.split)
    val <- ## mapply(fitfun,
        lapply(
            frm.split,fitfun,
            formula = as.formula(mform$model))

    use <- !sapply(val, is.null)
    if (nbad <- sum(!use))
        warning("Fitting failed for ",nbad," group(s), probably because a factor only had one level")

    ## Contrary to nlme, we keep the erronous ones as well
    pool <- !isGLM || .hasScale(family2char(family))
    new("lmList4", setNames(val, nms),
	call = mCall, pool = pool,
	groups = ordered(groups),
        origOrder = match(unique(as.character(groups)), nms)
        )
}

## (currently hidden) auxiliaries
isGLMlist <- function(object, ...) {
    D <- getDataPart(object)
    length(D) >= 1 && inherits(D[[1]], "glm")
}

## does a glm family have a "scale" [from stats:::logLik.glm() ] :
.hasScale <- function(family)
    family %in% c("gaussian", "Gamma", "inverse.gaussian")
family2char <- function(fam) {
    if(is.function(fam)) fam()$family else if(!is.character(fam)) fam$family else fam
}

##' Does a lmList4 object have a "scale" / sigma / useScale ?
hasScale <- function(object)
    !isGLMlist(object) || .hasScale(family(object[[1]])$family)


##' @importFrom stats coef
##' @S3method coef lmList4
## Extract the coefficients and form a  data.frame if possible
## FIXME: commented out nlme stuff (augFrame etc.).  Restore, or delete for good
## FIXME: modified so that non-estimated values will be NA rather than set to
##        coefs of first non-null estimate.  Is that OK??
coef.lmList4 <- function(object,
                        ## augFrame = FALSE, data = NULL,
                        ##which = NULL, FUN = mean, omitGroupingFactor = TRUE,
                        ...) {
    coefs <- lapply(object, coef)
    non.null <- !vapply(coefs, is.null, logical(1))
    if (any(non.null)) {
        template <- coefs[non.null][[1]]
        ## different parameter sets may be estimated for different subsets of data ...
        allnames <- Reduce(union, lapply(coefs[non.null], names))
        if (is.numeric(template)) {
            co <- matrix(NA,
                         ncol = length(allnames),
                         nrow = length(coefs),
                         dimnames = list(names(object), allnames))
            for (i in names(object)) {
                co[i,names(coefs[[i]])] <- coefs[[i]]
            }
            coefs <- as.data.frame(co)
            effectNames <- names(coefs)
            ## if(augFrame) {
            ##     if (is.null(data)) {
            ##         data <- getData(object)
            ##     }
            ##     data <- as.data.frame(data)
            ##     if (is.null(which)) {
            ##         which <- 1:ncol(data)
            ##     }
            ##     data <- data[, which, drop = FALSE]
            ##     ## eliminating columns with same names as effects
            ##     data <- data[, is.na(match(names(data), effectNames)), drop = FALSE]
            ##     data <- gsummary(data, FUN = FUN, groups = getGroups(object))
            ##     if (omitGroupingFactor) {
            ##         data <- data[, is.na(match(names(data),
            ##                                    names(getGroupsFormula(object,
            ##                                                           asList = TRUE)))),
            ##                      drop = FALSE]
            ##     }
            ##     if (length(data) > 0) {
            ##         coefs <- cbind(coefs, data[row.names(coefs),,drop = FALSE])
            ##     }
            ## }
            attr(coefs, "level") <- attr(object, "level")
            attr(coefs, "label") <- "Coefficients"
            attr(coefs, "effectNames") <- effectNames
            attr(coefs, "standardized") <- FALSE
        } ## is.numeric(template)
    }
    coefs
}

### FIXME?:  nlme *does* export this -- we export sigma()  [instead ?]
pooledSD <- function(x, allow.0.df = TRUE)
{
    stopifnot(is(x, "lmList4"))
    if(!hasScale(x)) {
        if(allow.0.df)
            return(structure(1, df = NA)) ## scale := 1  if(!useScale)
        ## else
        stop("no scale, hence no pooled SD for this object")
    }

    sumsqr <- rowSums(sapply(x,
                             function(el) {
                                 if (is.null(el)) {
                                     c(0,0)
                                 } else {
                                     res <- resid(el)
                                     c(sum(res^2), length(res) - length(coef(el)))
                                 }
			     }))
    if (sumsqr[2] == 0) { ## FIXME? rather return NA with a warning ??
        stop("No degrees of freedom for estimating std. dev.")
    }
    val <- sqrt(sumsqr[1]/sumsqr[2])
    attr(val, "df") <- sumsqr[2]
    val
}

sigma.lmList4 <- function(object, ...)
    if(hasScale(object)) as.vector(pooledSD(object)) else 1
## 1 for GLM  <==>  1 when useScale is FALSE for [G]LMMs


##' @importFrom methods show
##' @exportMethod show
setMethod("show", "lmList4", function(object)
{
    mCall <- object@call
    cat("Call:", deparse(mCall), "\n")
    cat("Coefficients:\n")
    print(coef(object))
    if (object@pool) {
        cat("\n")
        poolSD <- pooledSD(object)
        dfRes <- attr(poolSD, "df")
        RSE <- c(poolSD)
        cat("Degrees of freedom: ", length(unlist(lapply(object, fitted))),
            " total; ", dfRes, " residual\n", sep = "")
        cat("Residual standard error:", format(RSE))
        cat("\n")
    }
})

##' @S3method confint lmList4
confint.lmList4 <- function(object, parm, level = 0.95, ...)
{
    mCall <- match.call()
    if (length(object) < 1)
        return(new("lmList4.confint", array(numeric(0), c(0,0,0))))
    mCall$object <- object[[1]]
    ## the old recursive strategy doesn't work with S3 objects --
    ##  calls "confint.lmList4" again instead of calling "confint"
    mCall[[1]] <- quote(confint)
    template <- eval(mCall)
    if(is.null(d <- dim(template))) ## MASS:::confint.profile.glm() uses drop(), giving vector
	d <- dim(template <- rbind("(Intercept)" = template))
    template[] <- NA_real_
    val <- array(template, c(d, length(object)),
                 c(dimnames(template), list(names(object))))
    pool <- list(...)$pool
    if (is.null(pool)) pool <- object$pool
    if (length(pool) > 0 && pool[1]) { ## do our own
        sd <- pooledSD(object)
        a <- (1 - level)/2
        fac <- sd * qt(c(a, 1 - a)/2, attr(sd, "df"))
        parm <- dimnames(template)[[1]]
	for (i in seq_along(object))
	    if(!is.null(ob.i <- object[[i]]))
		val[ , , i] <- coef(ob.i)[parm] +
		    sqrt(diag(summary(object[[i]], corr = FALSE)$cov.unscaled
			      )[parm]) %o% fac
    } else { ## build on confint() method for "glm" / "lm" :
	for (i in seq_along(object))
	    if(!is.null(mCall$object <- object[[i]])) {
		ci <- eval(mCall)
		if(is.null(dim(ci))) ## MASS:::confint.profile.glm() ...
		    ci <- rbind("(Intercept)" = ci)
		if(identical(dim(ci), d))
		    val[ , , i] <- ci
		else ## some coefficients were not estimable
		    val[rownames(ci), , i] <- ci
	    }
    }
    new("lmList4.confint", aperm(val, 3:1))
}

##' @importFrom graphics plot
##' @importFrom lattice .......
##' @S3method plot lmList4.confint
plot.lmList4.confint <- function(x, y, order, ...)
{
##    stopifnot(require("lattice"))
    arr <- as(x, "array")
    dd <- dim(arr)
    dn <- dimnames(arr)
    levs <- dn[[1]]
    if (!missing(order) &&
        (ord <- round(order[1])) %in% seq(dd[3]))
        levs <- levs[order(rowSums(arr[ , , ord]))]
    ll <- length(arr)
    df <- data.frame(group =
                     ordered(rep(dn[[1]], dd[2] * dd[3]),
                             levels = levs),
                     intervals = as.vector(arr),
                     what = gl(dd[3], dd[1] * dd[2], length = ll, labels = dn[[3]]),
                     end = gl(dd[2], dd[1], length = ll))
    panelfun <- function(x, y, pch = dot.symbol$pch,
            col = dot.symbol$col, cex = dot.symbol$cex,
            font = dot.symbol$font, ...)
        {
            x <- as.numeric(x)
            y <- as.numeric(y)
            ok <- !is.na(x) & !is.na(y)
            yy <- y[ok]
            xx <- x[ok]
            dot.symbol <- trellis.par.get("dot.symbol")
            dot.line <- trellis.par.get("dot.line")
            panel.abline(h = yy, lwd = dot.line$lwd, lty = dot.line$lty, col =
                         dot.line$col)
            lpoints(xx, yy, pch = "|", col = col, cex = cex, font = font, ...)
            lower <- tapply(xx, yy, min)
            upper <- tapply(xx, yy, max)
            nams <- as.numeric(names(lower))
            lsegments(lower, nams, upper, nams, col = col, lty = 1, lwd =
                      if (dot.line$lwd) {
                          dot.line$lwd
                      } else {
                          2
                      })
        }
    dotplot(group ~ intervals | what,
            data = df,
            scales = list(x="free"),
            panel=panelfun, ...)
}

##' @importFrom stats update
##' @S3method update lmList4
update.lmList4 <- function(object, formula., ..., evaluate = TRUE) {
    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
        call$formula <- update.formula(formula(object), formula.)
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}

##' @importFrom stats formula
##' @S3method formula lmList4
##' @return of class "formula" ==> as.formula() rather than just [["formula"]]
formula.lmList4 <- function(x, ...) structure(x@call[["formula"]], class = "formula")

##' Get the grouping factor of an "lmList4" object
##' Important as auxiliary method for many of the nlme-imported methods:
getGroups.lmList4 <- function(object, ...) object@groups


### All the other "lmList4" S3 methods are imported from  nmle :
##
.ns.nlme <- asNamespace("nlme")
.ns.lme4 <- environment() ## == asNamespace("lme4") during build/load
##
## To do this, we need to make them use *our* namespace, e.g. to use our  pooledSD()
## However, then we get from codetools :
##
## fitted.lmList4: no visible global function definition for 'getGroups'
## pairs.lmList4: no visible global function definition for 'gsummary'
## pairs.lmList4: no visible global function definition for 'getGroups'
## plot.lmList4: no visible global function definition for 'c_deparse'
## plot.lmList4: no visible global function definition for 'getGroups'
## predict.lmList4: no visible global function definition for 'getGroups'
## print.lmList4: no visible global function definition for 'c_deparse'
## qqnorm.lmList4: no visible global function definition for 'getGroups'
## qqnorm.lmList4: no visible global function definition for 'gsummary'
## residuals.lmList4: no visible global function definition for 'getGroups'
##
## which we avoid via
for(fn in c("gsummary", "c_deparse")) {
    assign(fn, get(fn, envir = .ns.nlme, inherits=FALSE))
}

for(fn in c("fitted", "fixef", "logLik", "pairs", "plot", "predict",
            ## "print", <- have our own show()
           "qqnorm", "ranef", "residuals", "summary")) {
    meth <- get(paste(fn, "lmList",  sep="."), envir = .ns.nlme, inherits=FALSE)
    environment(meth) <- .ns.lme4 # e.g. in order to use *our* pooledSD()
    assign(paste(fn, "lmList4", sep="."), meth)
}
rm(fn)

