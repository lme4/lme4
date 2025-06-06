## copied/modified from nlme

##' split, on the nm call, the rhs of a formula into a list of subformulas
splitFormula <- function(form, sep = "/")
{
    if (inherits(form, "formula") ||
        mode(form) == "call" && form[[1]] == as.name("~"))
        splitFormula(form[[length(form)]], sep = sep)
    else if (mode(form) == "call" && form[[1]] == as.name(sep))
        do.call(c, lapply(as.list(form[-1]), splitFormula, sep = sep))
    else if (mode(form) == "(")
        splitFormula(form[[2]], sep = sep)
    else if (length(form))
        list(asOneSidedFormula(form))
    ## else
    ##  NULL
}

## Recursive version of all.vars
allVarsRec <- function(object)
{
    if (is.list(object)) {
        unlist(lapply(object, allVarsRec))
    } else {
        all.vars(object)
    }
}

## simple version of getData.gnls from nlme
## but we *should* and *can* work with environment(formula(.))
getData.merMod <-  function(object) {
    mCall <- getCall(object)
    data <- eval(mCall$data, environment(formula(object)))
    if (!is.data.frame(data) && !is.matrix(data)) stop(paste(sQuote("data"),"object found is not a data frame or matrix"))
    return(data)
}

asOneFormula <-
    ## Constructs a linear formula with all the variables used in a
    ## list of formulas, except for the names in omit
    function(..., omit = c(".", "pi"))
{
    names <- unique(allVarsRec(list(...)))
    names <- names[is.na(match(names, omit))]
    if (length(names))
        as.formula(paste("~", paste(names, collapse = "+"))) # else NULL
}

getIDLabels <- function(object, form=formula(object)) {
    mf <- factorize(form,model.frame(object))
    if (length(ff <- reformulas::findbars(form))>0) {
        grps <- lapply(ff,"[[",3)
    } else {
        grps <- form[[2]]
    }
    if (identical(grps,quote(.obs))) return(seq(fitted(object)))
    fList <- lapply(grps,function(x) eval(x,mf))
    do.call(interaction,fList)
}

## TESTING
## lme4:::getIDLabels(fm1)

## Return the formula(s) for the groups associated with object.
## The result is a one-sided formula unless asList is TRUE in which case
## it is a list of formulas, one for each level.
getGroupsFormula <- function(object, asList = FALSE, sep = "+")
    UseMethod("getGroupsFormula")

getGroupsFormula.default <-
    ## Return the formula(s) for the groups associated with object.
    ## The result is a one-sided formula unless asList is TRUE in which case
    ## it is a list of formulas, one for each level.
    function(object, asList = FALSE, sep = "/")
{
    form <- formula(object)
    if (!inherits(form, "formula")){
        stop("\"Form\" argument must be a formula")
    }
    form <- form[[length(form)]]
    if (!((length(form) == 3) && (form[[1]] == as.name("|")))) {
        ## no conditioning expression
        return(NULL)
    }
    ## val <- list( asOneSidedFormula( form[[ 3 ]] ) )
    val <- splitFormula(asOneSidedFormula(form[[3]]), sep = sep)
    names(val) <- unlist(lapply(val, function(el) deparse(el[[2]])))
                                        #  if (!missing(level)) {
                                        #    if (length(level) == 1) {
                                        #      return(val[[level]])
                                        #    } else {
                                        #      val <- val[level]
                                        #    }
                                        #  }
    if (asList) as.list(val)
    else as.formula(paste("~", paste(names(val), collapse = sep)))
}

getGroupsFormula.merMod <- function(object,asList=FALSE, sep="+") {
    if (asList) {
        lapply(names(object@flist),asOneSidedFormula)
    } else {
        asOneSidedFormula(paste(names(object@flist),collapse=sep))
    }
}

getCovariateFormula <- function (object)
{
    form <- formula(object)
    if (!(inherits(form, "formula"))) {
        stop("formula(object) must return a formula")
    }
    form <- form[[length(form)]]
    if (length(form) == 3 && form[[1]] == as.name("|")) {
        form <- form[[2]]
    }
    eval(substitute(~form))
}

getResponseFormula <-
    function(object)
{
    ## Return the response formula as a one sided formula
    form <- formula(object)
    if (!(inherits(form, "formula") && (length(form) == 3))) {
        stop("\"Form\" must be a two sided formula")
    }
    as.formula(paste("~", deparse(form[[2]])))
}

##' diagnostic plots for merMod fits
##' @param x a fitted [ng]lmer model
##' @param form an optional formula specifying the desired type of plot. Any
##' variable present in the original data frame used to obtain
##' \code{x} can be referenced. In addition, \code{x} itself can be
##' referenced in the formula using the symbol \code{"."}. Conditional
##'  expressions on the right of a \code{|} operator can be used to
##'  define separate panels in a lattice display. Default is
##' \code{resid(., type = "pearson") ~ fitted(.)}, corresponding to a plot
##' of the standardized residuals versus fitted values.
##' @param abline an optional numeric value, or numeric vector of length
##'   two. If given as a single value, a horizontal line will be added to the
##'   plot at that coordinate; else, if given as a vector, its values are
##'   used as the intercept and slope for a line added to the plot. If
##'   missing, no lines are added to the plot.
##' @param id an optional numeric value, or one-sided formula. If given as
##' a value, it is used as a significance level for a two-sided outlier
##' test for the standardized, or normalized residuals. Observations with
##'   absolute standardized (normalized) residuals greater than the \eqn{1-value/2}
##' quantile of the standard normal distribution are
##' identified in the plot using \code{idLabels}. If given as a one-sided
##'   formula, its right hand side must evaluate to a  logical, integer, or
##'   character vector which is used to identify observations in the
##'   plot. If missing, no observations are identified.
##' @param idLabels an optional vector, or one-sided formula. If given as a
##'   vector, it is converted to character and used to label the
##'   observations identified according to \code{id}. If given as a
##'    vector, it is converted to character and used to label the
##'    observations identified according to \code{id}. If given as a
##'    one-sided formula, its right hand side must evaluate to a vector
##'    which is converted to character and used to label the identified
##'    observations. Default is the interaction of all the grouping variables
##'    in the data frame.  The special formula
##' @param grid an optional logical value indicating whether a grid should
##'    be added to plot. Default depends on the type of lattice plot used:
##'    if \code{xyplot} defaults to \code{TRUE}, else defaults to
##'    \code{FALSE}.
##'  @param \dots optional arguments passed to the lattice plot function.
##' @details Diagnostic plots for the linear mixed-effects fit are obtained. The
##'  \code{form} argument gives considerable flexibility in the type of
##'  plot specification. A conditioning expression (on the right side of a
##'  \code{|} operator) always implies that different panels are used for
##'  each level of the conditioning factor, according to a lattice
##'  display. If \code{form} is a one-sided formula, histograms of the
##'  variable on the right hand side of the formula, before a \code{|}
##'  operator, are displayed (the lattice function \code{histogram} is
##'  used). If \code{form} is two-sided and both its left and
##'  right hand side variables are numeric, scatter plots are displayed
##'  (the lattice function \code{xyplot} is used). Finally, if \code{form}
##'  is two-sided and its left had side variable is a factor, box-plots of
##'  the right hand side variable by the levels of the left hand side
##'  variable are displayed (the lattice function  \code{bwplot} is used).
##' @author original version in \code{nlme} package by Jose Pinheiro and Douglas Bates
##' @examples
##' data(Orthodont,package="nlme")
##' fm1 <- lmer(distance ~ age + (age|Subject), data=Orthodont)
##' ## standardized residuals versus fitted values by gender
##' plot(fm1, resid(., scaled=TRUE) ~ fitted(.) | Sex, abline = 0)
##' ## box-plots of residuals by Subject
##' plot(fm1, Subject ~ resid(., scaled=TRUE))
##' ## observed versus fitted values by Subject
##' plot(fm1, distance ~ fitted(.) | Subject, abline = c(0,1))
##' ## residuals by age, separated by Subject
##' plot(fm1, resid(., scaled=TRUE) ~ age | Sex, abline = 0)

##' if (require(ggplot2)) {
##'     ## we can create the same plots using ggplot2 and the fortify() function
##'     fm1F <- fortify(fm1)
##'     ggplot(fm1F, aes(.fitted,.resid)) + geom_point(colour="blue") +
##'            facet_grid(.~Sex) + geom_hline(yintercept=0)
##'     ## note: Subjects are ordered by mean distance
##'     ggplot(fm1F, aes(Subject,.resid)) + geom_boxplot() + coord_flip()
##'     ggplot(fm1F, aes(.fitted,distance))+ geom_point(colour="blue") +
##'         facet_wrap(~Subject) +geom_abline(intercept=0,slope=1)
##'     ggplot(fm1F, aes(age,.resid)) + geom_point(colour="blue") + facet_grid(.~Sex) +
##'         geom_hline(yintercept=0)+geom_line(aes(group=Subject),alpha=0.4)+geom_smooth(method="loess")
##'     ## (warnings about loess are due to having only 4 unique x values)
##'     detach("package:ggplot2")
##' }
##' @S3method plot merMod
##' @method plot merMod
##' @export
plot.merMod <-
    function(x, form = resid(., type = "pearson") ~ fitted(.), abline,
             id = NULL, idLabels = NULL,
             grid, ...)
        ## Diagnostic plots based on residuals and/or fitted values
{
    object <- x
    if (!inherits(form, "formula"))
        stop("\"form\" must be a formula")
    ## constructing data
    ## can I get away with using object@frame???
    allV <- all.vars(asOneFormula(form, id, idLabels))
    allV <- allV[is.na(match(allV,c("T","F","TRUE","FALSE",".obs")))]
    if (length(allV) > 0) {
        data <- getData(object)
        if (is.null(data)) {            # try to construct data
            alist <- lapply(as.list(allV), as.name)
            names(alist) <- allV
            alist <- c(list(as.name("data.frame")), alist)
            mode(alist) <- "call"
            data <- eval(alist, sys.parent(1))
        } else if (any(naV <- is.na(match(allV, names(data)))))
              stop(allV[naV], " not found in data")
    } else data <- NULL

    ## this won't do because there may well be variables we want
    ##  that were not in the model call

    ## data <- object@frame

    ## argument list
    dots <- list(...)
    args <- if (length(dots) > 0) dots else list()
    ## appending object to data, and adding observation-number variable
    if (length(data) > 0) {
        data <- cbind(data, .obs = seq(nrow(data)))
    }
    data <- as.list(c(as.list(data), . = list(object)))
    ## covariate - must always be present
    covF <- getCovariateFormula(form)
    .x <- eval(covF[[2]], data)
    if (!is.numeric(.x)) {
        stop("Covariate must be numeric")
    }
    argForm <- ~ .x
    argData <- data.frame(.x = .x, check.names = FALSE)
    if (is.null(args$xlab)) {
        if (is.null(xlab <- attr(.x, "label")))
            xlab <- deparse(covF[[2]])
        args$xlab <- xlab
    }

    ## response - need not be present
    respF <- getResponseFormula(form)
    if (!is.null(respF)) {
        .y <- eval(respF[[2]], data)
        if (is.null(args$ylab)) {
            if (is.null(ylab <- attr(.y, "label")))
                ylab <- deparse(respF[[2]])
            args$ylab <- ylab
        }
        argForm <- .y ~ .x
        argData[, ".y"] <- .y
    }

    ## groups - need not be present
    grpsF <- getGroupsFormula(form)
    if (!is.null(grpsF)) {
        ## ?? FIXME ???
        gr <- splitFormula(grpsF, sep = "*")
        for(i in seq_along(gr)) {
            auxGr <- all.vars(gr[[i]])
            for(j in auxGr)
                argData[[j]] <- eval(as.name(j), data)
        }
        argForm <-
            as.formula(paste(if (length(argForm) == 2)
                                 "~ .x |" else ".y ~ .x |",
                             deparse(grpsF[[2]])))
    }
    ## adding to args list
    args <- c(list(argForm, data = argData), args)
    if (is.null(args$strip)) {
        args$strip <- function(...) strip.default(..., style = 1)
    }
    if (is.null(args$cex)) args$cex <- par("cex")
    if (is.null(args$adj)) args$adj <- par("adj")

    if (!is.null(id)) {       ## identify points in plot
        idResType <- "pearson"  ## diff from plot.lme: 'normalized' not available
        id <- switch(mode(id),
                     numeric = {
                         if (id <= 0 || id >= 1)
                             stop(shQuote("id")," must be between 0 and 1")
                         abs(resid(object, type = idResType))/sigma(object) >
                             -qnorm(id / 2)
                     },
                     call = eval(asOneSidedFormula(id)[[2]], data),
                     stop(shQuote("id")," can only be a formula or numeric.")
                     )
        if (is.null(idLabels)) {
            idLabels <- getIDLabels(object)
        } else {
            if (inherits(idLabels,"formula")) {
                idLabels <- getIDLabels(object,idLabels)
            } else if (is.vector(idLabels)) {
                if (length(idLabels <- unlist(idLabels)) != length(id)) {
                    stop("\"idLabels\" of incorrect length")
                }
            } else stop("\"idLabels\" can only be a formula or a vector")
        }
        ## DON'T subscript by id, will be done later
        idLabels <- as.character(idLabels)
    }

    ## defining abline, if needed
    if (missing(abline)) {
        abline <- if (missing(form)) # r ~ f
                      c(0, 0) else NULL
    }

                                        #assign("id", id , where = 1)
                                        #assign("idLabels", idLabels, where = 1)
                                        #assign("abl", abline, where = 1)
    assign("abl", abline)

    ## defining the type of plot
    if (length(argForm) == 3) {
        if (is.numeric(.y)) {           # xyplot
            plotFun <- "xyplot"
            if (is.null(args$panel)) {
                args <- c(args,
                          panel = list(function(x, y, subscripts, ...)
                          {
                              x <- as.numeric(x)
                              y <- as.numeric(y)
                              dots <- list(...)
                              if (grid) panel.grid()
                              panel.xyplot(x, y, ...)
                              if (any(ids <- id[subscripts])){
                                  ltext(x[ids], y[ids], idLabels[subscripts][ids],
                                        cex = dots$cex, adj = dots$adj)
                              }
                              if (!is.null(abl)) {
                                  if (length(abl) == 2) panel.abline(a = abl, ...) else panel.abline(h = abl, ...)
                              }
                          }))
            }
        } else {                                # assume factor or character
            plotFun <- "bwplot"
            if (is.null(args$panel)) {
                args <- c(args,
                          panel = list(function(x, y, ...)
                          {
                              if (grid) panel.grid()
                              panel.bwplot(x, y, ...)
                              if (!is.null(abl)) {
                                  panel.abline(v = abl[1], ...)
                              }
                          }))
            }
        }
    } else {
        plotFun <- "histogram"
        if (is.null(args$panel)) {
            args <- c(args,
                      panel = list(function(x, ...)
                      {
                          if (grid) panel.grid()
                          panel.histogram(x, ...)
                          if (!is.null(abl)) {
                              panel.abline(v = abl[1], ...)
                          }
                      }))
        }
    }

    ## defining grid
    if (missing(grid)) {
        grid <- (plotFun == "xyplot")
    }
                                        # assign("grid", grid, where = 1)
    do.call(plotFun, as.list(args))
}

## no longer defining `fortify` S3 generic

##' @rdname fortify
##' @S3method fortify lmerMod
##' @method fortify lmerMod
##' @export
##'   as function, not as S3 method, see ../man/fortify.Rd  :
fortify.merMod <- function(model, data=getData(model), ...) {

    ## FIXME: get influence measures via influence.ME?
    ##   (expensive, induces dependency ...)
    ## FIXME: different kinds of residuals?
    ## FIXME: deal with na.omit/predict etc.
    data$.fitted <- predict(model)
    data$.resid <- resid(model)
    data$.scresid <- resid(model,type="pearson",scaled=TRUE)
    data
}


## autoplot???

##  plot method for plot.summary.mer ... coefplot-style
##  horizontal, vertical? other options???
##  scale?
plot.summary.mer <- function(object, type="fixef", ...) {
    if(any(!type %in% c("fixef","vcov")))
        stop("'type' not yet implemented: ", type)
    stop("FIXME -- not yet implemented")

}

## TO DO: allow faceting formula
## TO DO: allow qqline to be optional
## TO DO (harder): steal machinery from qq.gam for better GLMM Q-Q plots
qqmath.merMod <- function(x, data = NULL, id=NULL, idLabels=NULL, ...) {
    ## klugey attempt to detect whether user forgot to specify argument
    ## names explicitly (after addition of required 'data' argument)
    ## NOT completely tested!
    if (!is.null(data)) {
        idLabels <- id
        id <- data
        warning("qqmath.merMod takes ", sQuote("data"), "as its ",
                "first argument for S3 method compatibility: ",
                "in the future, please ",
                "specify the ", sQuote("id"), " and ",
                sQuote("idLabels"), " arguments explicitly ",
                "i.e. ",
                sQuote("qqmath(fitted_model, id = ..., [idLabels = ...], ...)"))
    }
   
    values <- residuals(x, type="pearson", scaled=TRUE)
    data <- getData(x)
    ## DRY: copied from plot.merMod, should modularize/refactor
    if (!is.null(id)) {       ## identify points in plot
        id <- switch(mode(id),
                     numeric = {
                         if (id <= 0 || id >= 1)
                             stop(shQuote("id")," must be between 0 and 1")
                         as.logical(abs(values) > -qnorm(id / 2))
                     },
                     call = eval(asOneSidedFormula(id)[[2]], data),
                     stop(shQuote("id")," can only be a formula or numeric.")
                     )
        if (is.null(idLabels)) {
            idLabels <- getIDLabels(x)
        } else {
            if (inherits(idLabels,"formula")) {
                idLabels <- getIDLabels(x,idLabels)
            } else if (is.vector(idLabels)) {
                if (length(idLabels <- unlist(idLabels)) != length(id)) {
                    stop("\"idLabels\" of incorrect length")
                }
            } else stop("\"idLabels\" can only be a formula or a vector")
        }
        idLabels <- as.character(idLabels)

    }
    ## DON'T subscript by id, will be done later
    qqpanel <- function(x, subscripts, ...) {
        dots <- list(...)
        panel.qqmathline(x, ...)
        panel.qqmath(x, ...)
        if (any(ids <- id[subscripts])) {
            xs <- x[subscripts]
            pp <- setNames(ppoints(length(xs)),
                           names(sort(xs)))
            ## want to plot qnorm(pp) vs sort(x)
            ## ... but want to pick out the elements that corresponded
            ## to ids **before** sorting
            xx <- qnorm(pp)[names(xs)[ids]]
            yy <- sort(x)[names(xs)][ids] ## quantile(values, pp)[ids]
            ltext(xx,
                  yy,
                  idLabels[ids],
                  cex = dots$cex, adj = dots$adj)
        }
    }
    qqmath(values, xlab = "Standard normal quantiles",
           ylab = "Standardized residuals",
           prepanel = prepanel.qqmathline,
           panel = qqpanel,
           ...)
}

## qqmath(~residuals(gm1)|cbpp$herd)
