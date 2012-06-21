##' Fit a linear mixed model (LMM)
##'
##' @title Fit Linear Mixed-Effects Models
##' @concept LMM
##' @aliases lmer
##' @param formula a two-sided linear formula object describing the
##'    fixed-effects part of the model, with the response on the left of a
##'    \code{~} operator and the terms, separated by \code{+} operators, on
##'    the right.  The vertical bar character \code{"|"} separates an
##'    expression for a model matrix and a grouping factor.
##' @param data an optional data frame containing the variables named in
##'    \code{formula}.  By default the variables are taken from the environment
##'    from which \code{lmer} is called.
##' @param REML logical scalar - Should the estimates be chosen to optimize
##'    the REML criterion (as opposed to the log-likelihood)?  Defaults to
##'    \code{TRUE}.
##' @param sparseX logical - should a sparse model matrix be used for the
##'    fixed-effects terms?  Defaults to \code{FALSE}. Currently inactive.
##' @param control a named list of control parameters for the estimation
##'     algorithm, specifying only the ones to be changed from their
##'     default values.  Hence defaults to an empty list.\cr
##'     Possible control options and their default values are:
##'   \describe{
##'      \item{\code{msVerbose}:}{a logical value passed as the
##'      \code{trace} argument to \code{nlminb} (see documentation on
##'      that function).  Default is \code{getOption("verbose")}.}
##'	\item{\code{maxIter}:}{a positive integer passed as the
##'	\code{maxIter} argument to \code{nlminb} (see documentation on
##'	that function).  Default is \code{300}.}
##'	\item{\code{maxFN}:}{a positive integer specifying the
##'	 maximum number of evaluations of the deviance function allowed
##'	 during the optimization. Default is \code{900}.}
##'	\item{\code{tol}:}{a positive number specifying the
##'	 convergence tolerance, currently only for the PWRSS iterations
##'	 in \code{\link{glmer}}.  Default is \code{0.000001}.}
##'   }
##' @param start a named list of starting values for the parameters in the
##'    model.  For \code{lmer} this can be a numeric vector or a list with one
##'    component named \code{"theta"}. Infrequently used.
##' @param verbose integer scalar.  If \code{> 0} verbose output is generated
##'    during the optimization of the parameter estimates.  If \code{> 1} verbose
##'    output is generated during the individual PIRLS steps.
##' @param subset an optional expression indicating the subset of the rows of
##'     \code{data} that should be used in the fit. This can be a logical
##'     vector, or a numeric vector indicating which observation numbers are
##'     to be included, or a  character  vector of the row names to be
##'     included.  All observations are included by default.
##' @param weights an optional vector of \sQuote{prior weights} to be used in the
##'     fitting process.  Should be \code{NULL} or a numeric vector.
##' @param na.action a function that indicates what should happen when the
##'     data contain \code{NA}s.  The default action (\code{na.fail}) prints
##'     an error message and terminates if there are any incomplete
##'     observations.
##' @param offset this can be used to specify an \emph{a priori} known component
##'     to be included in the linear predictor during fitting. This should be
##'     \code{NULL} or a numeric vector of length equal to the number of cases.
##'     One or more \code{\link{offset}} terms can be included in the formula
##'     instead or as well, and if more than one is specified their sum is used.
##'     See \code{\link{model.offset}}.
##' @param contrasts an optional list. See the \code{contrasts.arg} of
##'     \code{model.matrix.default}.
##' @param devFunOnly logical - return only the deviance evaluation function.
##' @param optimizer character - name of optimizing function
##' @param \dots other potential arguments.  A \code{method} argument was used
##'    in earlier versions of the package. Its functionality has been replaced by
##'    the \code{REML} argument.
##' @return An object of class \code{"\linkS4class{merMod}"}, for which many
##'    methods are available.  See there for details.
##' @seealso The \code{\linkS4class{merMod}} class, \code{\link[stats]{lm}}
##' @keywords models
##' @examples
##' ## linear mixed models - reference values from older code
##' (fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
##' (fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' anova(fm1, fm2)
##' @export
##' @importFrom minqa bobyqa
lmer <- function(formula, data, REML = TRUE, sparseX = FALSE,
                 control = list(), start = NULL,
                 verbose = 0L, subset, weights, na.action, offset,
                 contrasts = NULL, devFunOnly=FALSE,
                 optimizer="Nelder_Mead", ...)
{
    if (sparseX) warning("sparseX = TRUE has no effect at present")
    mf <- mc <- match.call()
    ## '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
        if (!is.null(l...$family)) {  # call glmer if family specified
            warning("calling lmer with family() is deprecated: please use glmer() instead")
            mc[[1]] <- as.name("glmer")
            return(eval(mc, parent.frame()) )
        }
        ## Check for method argument which is no longer used
        if (!is.null(method <- l...$method)) {
            msg <- paste("Argument", sQuote("method"), "is deprecated.")
            ## FIXME: this will fail if method *not* in ("Laplace","AGQ") ...
            if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
                warning(msg)
                l... <- l...[names(l...) != "method"]
            } else stop(msg)
        }
        if(length(l...))
            warning("extra argument(s) ",
                    paste(sQuote(names(l...)), collapse=", "),
                    " disregarded")
    }

    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
                                        # evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substituted "|" by "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
                                        # random effects and terms modules
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    if (FALSE) {
        ## BMB: I clearly don't know what's going on here yet.
        ## test: lmer(angle ~ temp + recipe + (1 | replicate), data = cake)
        ## test: sstudy9 <- subset(sleepstudy, Days == 1 | Days == 9)
        ## m1 <- lmer(Reaction ~ 1 + Days + (1 + Days | Subject), data = sstudy9)
        rankZ1 <- rankMatrix(bdiag(reTrms$Zt,Diagonal(ncol(reTrms$Zt))))
        pZ1 <- nrow(reTrms$Zt)+ncol(reTrms$Zt)
        if (rankZ1<pZ1)
            stop(gettextf("rank of cBind(Z,1) = %d < ncol(Z)+1 = %d", rankZ1, pZ1+1))
    }
    if (any(unlist(lapply(reTrms$flist, nlevels)) >= nrow(fr)))
        stop("number of levels of each grouping factor must be ",
             "less than number of obs")
    ## fixed-effects model matrix X - remove random effects from formula:
    form <- formula
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.matrix(form, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    if ((qrX <- qr(X))$rank < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    rho <- new.env(parent=parent.env(environment()))
    rho$pp <- do.call(merPredD$new, c(reTrms[c("Zt","theta","Lambdat","Lind")], n=nrow(X), list(X=X)))
    rho$resp <- mkRespMod(fr, if(REML) p else 0L)

    devfun <- mkdevfun(rho, 0L)
    devfun(reTrms$theta) # one evaluation to ensure all values are set

    if (devFunOnly) return(devfun)

    opt <- optwrap(optimizer,
                   devfun, rho$pp$theta, lower=reTrms$lower, control=control,
                   adj=FALSE)

    mkMerMod(environment(devfun), opt, reTrms, fr, mc)
}## { lmer }


##' Fit a generalized linear mixed model (GLMM)
##'
##' Fit a generalized linear mixed model, which incorporates both fixed-effects
##' parameters and random effects in a linear predictor, via maximum
##' likelihood. The linear predictor is related to the conditional
##' mean of the response through the inverse link function defined in
##' the GLM \code{family}.
##'
##' The expression for the likelihood of a mixed-effects model is an integral over
##' the random effects space. For a linear mixed-effects model (LMM), as fit by
##' \code{\link{lmer}}, this integral can be evaluated exactly.  For a
##' GLMM the integral must be approximated.  The most reliable
##' approximation for GLMMs with a single grouping factor for the
##' random effects is adaptive Gauss-Hermite quadrature.  The
##' \code{nAGQ} argument controls the number of nodes in the
##' quadrature formula.  A model with a single, scalar random-effects
##' term could reasonably use up to 25 quadrature points per scalar
##' integral.
##'
##' With vector-valued random effects the complexity of the
##' Gauss-Hermite quadrature formulas increases dramatically with the
##' dimension.  For a 3-dimensional vector-valued random effect
##' \code{nAGQ=5} requires 93 evaluations of the GLM deviance per
##' evaluation of the approximate GLMM deviance.  For 20-dimensional
##  vector-valued random effects, \code{nAGQ=2} requires 41
##' evaluations of the GLM deviance per evaluation of the approximate
##' GLMM deviance.
##'
##' The default approximation is the Laplace approximation,
##' corresponding to \code{nAGQ=1}.
##' @title Fit Generalized Linear Mixed-Effects Models
##' @concept GLMM
##' @param family a GLM family, see \code{\link[stats]{glm}} and
##'    \code{\link[stats]{family}}.
##' @param compDev logical scalar - should compiled code be used for the
##'    deviance evaluation during the optimization of the parameter
##'    estimates?  Defaults to \code{TRUE}.
##' @param nAGQ integer scalar - the number of points per axis for evaluating
##'    the adaptive Gauss-Hermite approximation to the log-likelihood.  Applies
##'    only to \code{glmer} and defaults to 1, corresponding to the Laplace
##'    approximation.  Values greater than 1 produce greater accuracy in
##'    the evaluation of the log-likelihood at the expense of speed.  A value
##'    of zero use a faster but less exact form of parameter estimation for
##'    GLMMs by optimizing the random effects and the fixed-effects coefficients
##'    in the penalized iteratively reweighted least squares step.
##' @param start a named list of starting values for the parameters in the
##'    model.  If the list contains components named \code{fixef} and/or
##'    \code{theta}, these are used as the starting values for those slots.
##'    A numeric \code{start} argument of the appropriate length is used as the
##'    starting value of \code{theta}.
##' @param optimizer which optimizer(s) to use for each phase of optimization.
##'    A character vector or list of functions.
##'    If \code{length(optimizer)==2}, the first element will be used
##'    for the preliminary (random effects parameters only) optimization, while
##'    the second will be used for the final (random effects plus
##'    fixed effect parameters) phase. The built-in optimizers are
##'    \code{\link{Nelder_Mead}} and \code{\link[minqa]{bobyqa}} (from
##'    the \pkg{minqa} package; the default
##'    is to use \code{\link[minqa]{bobyqa}} for the first and
##'    \code{\link{Nelder_Mead}} for the final phase.
##'    (FIXME: simplify if possible!). For difficult model fits we have found
##'    \code{\link{Nelder_Mead}} to be more reliable but occasionally slower than
##'    \code{\link{bobyqa}}. Any minimizing function that allows
##'    box constraints can be used
##'    provided that it (1) takes input parameters \code{fn} (function
##'    to be optimized), \code{par} (starting parameter values),
##'    \code{lower} (lower bounds) and \code{control} (control parameters,
##'    passed through from the \code{control} argument) and (2)
##'    returns a list with (at least) elements \code{par} (best-fit
##'    parameters), \code{fval} (best-fit function value), \code{conv}
##'    (convergence code) and (optionally) \code{message} (informational
##'    message, or explanation of convergence failure).
##'    Special provisions are made for \code{\link{bobyqa}},
##'    \code{\link{Nelder_Mead}}, and optimizers wrapped in
##'    the \pkg{optimx} package; to use \pkg{optimx} optimizers
##'    (including \code{L-BFGS-B} from base \code{\link{optim}} and
##'    \code{\link{nlminb}}), pass the \code{method} argument to \code{optim}
##'    in the \code{control} argument.
##'
##' @param mustart optional starting values on the scale of the conditional mean,
##'    as in \code{\link[stats]{glm}}; see there for details.
##' @param etastart optional starting values on the scale of the unbounded
##'    predictor as in \code{\link[stats]{glm}}; see there for details.
##' @param tolPwrss numeric scalar - the tolerance for declaring convergence in
##'    the penalized iteratively weighted residual sum-of-squares step.
##'    Defaults to 1e-10.
##' @param \dots other potential arguments.  A \code{method} argument was used
##'    in earlier versions of the package. Its functionality has been replaced by
##'    the \code{nAGQ} argument.
##' @inheritParams lmer
##' @return An object of class \code{"\linkS4class{merMod}"}, for which many
##'    methods are available.  See there for details.
##' @seealso The \code{\linkS4class{merMod}} class, \code{\link[stats]{glm}}
##' @keywords models
##' @examples
##' ## generalized linear mixed model
##' library(lattice)
##' xyplot(incidence/size ~ period|herd, cbpp, type=c('g','p','l'),
##'        layout=c(3,5), index.cond = function(x,y)max(y))
##' (gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'               data = cbpp, family = binomial))
##' ## using nAGQ=0 only gets close to the optimum
##' (gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'                cbpp, binomial, nAGQ = 0))
##' ## using  nAGQ = 9  provides a better evaluation of the deviance
##' ## Currently the internal calculations use the sum of deviance residuals,
##' ## which is not directly comparable with the nAGQ=0 or nAGQ=1 result.
##' (gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'                cbpp, binomial, nAGQ = 9))
##'
##' ## GLMM with individual-level variability (accounting for overdispersion)
##' ## For this data set the model is the same as one allowing for a period:herd
##' ## interaction, which the plot indicates could be needed.
##' cbpp$obs <- 1:nrow(cbpp)
##' (gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
##'     (1 | herd) +  (1|obs),
##'               family = binomial, data = cbpp))
##' anova(gm1,gm2)
##' @export
glmer <- function(formula, data, family = gaussian, sparseX = FALSE,
                  control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                  compDev = TRUE, subset, weights, na.action, offset,
                  contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
                  tolPwrss = 1e-7, optimizer=c("bobyqa","Nelder_Mead"), ...)
{
    verbose <- as.integer(verbose)
    mf <- mc <- match.call()
                                        # extract family, call lmer for gaussian
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if( is.function(family)) family <- family()
    if (isTRUE(all.equal(family, gaussian()))) {
        mc[[1]] <- as.name("lmer")
        mc["family"] <- NULL            # to avoid an infinite loop
        return(eval(mc, parent.frame()))
    }

    if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
        stop('"quasi" families cannot be used in glmer')
### '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
        ## Check for invalid specifications
        if (!is.null(method <- list(...)$method)) {
            msg <- paste("Argument", sQuote("method"),
                         "is deprecated.\nUse", sQuote("nAGQ"),
                         "to choose AGQ.  PQL is not available.")
            if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
                warning(msg)
                l... <- l...[names(l...) != "method"]
            } else stop(msg)
        }
        if(length(l...))
            warning("extra argument(s) ",
                    paste(sQuote(names(l...)), collapse=", "),
                    " disregarded")
    }

    stopifnot(length(nAGQ <- as.integer(nAGQ)) == 1L,
              nAGQ >= 0L,
              nAGQ <= 25L,
              length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
                                        # evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "mustart", "etastart"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" for "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
                                        # random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    if ((maxlevels <- max(unlist(lapply(reTrms$flist, nlevels)))) > nrow(fr))
        stop("number of levels of each grouping factor must be",
             "greater than or equal to number of obs")
    ## FIXME: adjust test for families with estimated scale;
    ##   useSc is not defined yet/not defined properly?
    ##  if (useSc && maxlevels == nrow(fr))
    ##          stop("number of levels of each grouping factor must be",
    ##                "greater than number of obs")

    ## fixed-effects model matrix X - remove random parts from formula:
    form <- formula
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.matrix(form, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)

    ## Environment for deviance function.  For the optimizers the
    ## deviance function must be a simple function of a numeric
    ## parameter.  We put all the other information in the
    ## environment rho which is assigned as the environment of the
    ## deviance function.
    rho             <- as.environment(list(verbose=verbose, tolPwrss=tolPwrss))
    parent.env(rho) <- parent.frame()
    rho$pp          <- do.call(merPredD$new,
                               c(reTrms[c("Zt","theta","Lambdat","Lind")],
                                 n=nrow(X), list(X=X)))
    rho$resp        <- mkRespMod(fr, family=family)
    if (length(unique(rho$resp$y)) < 2L)
        stop("Response is constant - cannot fit the model")
    rho$verbose     <- as.integer(verbose)
                                        # initialize (from mustart)
    .Call(glmerLaplace, rho$pp$ptr(), rho$resp$ptr(), 0L, tolPwrss, verbose)
    rho$lp0         <- rho$pp$linPred(1) # each pwrss opt begins at this eta
    rho$pwrssUpdate <- glmerPwrssUpdate
    rho$compDev     <- compDev
    rho$lower       <- reTrms$lower     # not needed in rho?
    devfun <- mkdevfun(rho, 0L)
    if (devFunOnly && !nAGQ) return(devfun)
    if (length(optimizer)==1) {
        optimizer <- replicate(2,optimizer)
    }
    opt <- optwrap(optimizer[[1]],devfun,rho$pp$theta, rho$lower,
                   control=control,
                   adj=FALSE, verbose=verbose)
    rho$control <- attr(opt,"control")

    rho$nAGQ <- nAGQ
    if (nAGQ > 0L) {
        rho$nAGQ       <- nAGQ
        rho$lower      <- c(rho$lower, rep.int(-Inf, length(rho$pp$beta0)))
        rho$lp0        <- rho$pp$linPred(1)
        rho$dpars      <- seq_along(rho$pp$theta)
        rho$baseOffset <- rho$resp$offset + 0 # forcing a copy (!)
        rho$GQmat      <- GHrule(nAGQ)
        rho$fac        <- reTrms$flist[[1]]
        if (nAGQ > 1L) {
            if (length(reTrms$flist) != 1L || length(reTrms$cnms[[1]]) != 1L)
                stop("nAGQ > 1 is only available for models with a single, scalar random-effects term")
        }
        devfun <- mkdevfun(rho, nAGQ)
        if (devFunOnly) return(devfun)

        opt <- optwrap(optimizer[[2]],devfun,c(rho$pp$theta, rho$pp$delb),
                       rho$lower, control=control,
                       adj=TRUE, verbose=verbose)
        rho$resp$setOffset(rho$baseOffset)
    }
    mkMerMod(environment(devfun), opt, reTrms, fr, mc)
}## {glmer}

##' Fit a nonlinear mixed-effects model
##'
##' Fit nonlinear mixed-effects models, such as those used in
##' population pharmacokinetics.
##' @title Fit Nonlinear Mixed-Effects Models
##' @param formula a nonlinear mixed model formula (see detailed documentation)
##' @param start starting estimates for the nonlinear model
##'    parameters, as a named numeric vector or as a list with components
##'    \describe{
##'    \item{nlpars}{required numeric vector of starting values for the
##'         nonlinear model parameters}
##'    \item{theta}{optional numeric vector of starting values for the
##'         covariance parameters}
##'    }
##' @param \dots other potential arguments.  A \code{method} argument was used
##'    in earlier versions of the package. Its functionality has been replaced by
##'    the \code{nAGQ} argument.
##' @inheritParams glmer
##' @keywords models
##' @examples
##' ## nonlinear mixed models --- 3-part formulas ---
##'
##' (nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
##'              Orange, start = c(Asym = 200, xmid = 725, scal = 350)))
##' (nm1a <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
##'               Orange, start = c(Asym = 200, xmid = 725, scal = 350),
##'               nAGQ = 0L))
##' @export
nlmer <- function(formula, data, control = list(), start = NULL, verbose = 0L,
                  nAGQ = 1L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly = 0L, tolPwrss = 1e-10,
                  optimizer="Nelder_Mead", ...)
{
    vals <- nlformula(mc <- match.call())
    if ((qrX <- qr(X <- vals$X))$rank < (p <- ncol(X)))
        stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    rho <- list2env(list(verbose=verbose,
                         tolPwrss=0.001, # this is reset to the tolPwrss argument's value later
                         resp=vals$resp,
                         lower=vals$reTrms$lower),
                    parent=parent.frame())
    rho$pp <- do.call(merPredD$new,
                      c(vals$reTrms[c("Zt","theta","Lambdat","Lind")],
                        list(X=X, n=length(vals$respMod$mu), Xwts=vals$respMod$sqrtXwt,
                             beta0=qr.coef(qrX, unlist(lapply(vals$pnames, get,
                             envir = rho$resp$nlenv))))))
    rho$u0 <- rho$pp$u0
    rho$beta0 <- rho$pp$beta0
    devfun <- mkdevfun(rho, 0L) # deviance as a function of theta only
    if (devFunOnly && !nAGQ) return(devfun)
    devfun(rho$pp$theta) # initial coarse evaluation to get u0 and beta0
    rho$u0 <- rho$pp$u0
    rho$beta0 <- rho$pp$beta0
    rho$tolPwrss <- tolPwrss # Resetting this is intentional. The initial optimization is coarse.

    if (length(optimizer)==1) {
        optimizer <- replicate(2,optimizer)
    }

    opt <- optwrap(optimizer[[1]], devfun, rho$pp$theta, rho$lower,
                   control=control, adj=FALSE)
    rho$control <- attr(opt,"control")

    if (nAGQ > 0L) {
        rho$lower <- c(rho$lower, rep.int(-Inf, length(rho$beta0)))
        rho$u0    <- rho$pp$u0
        rho$beta0 <- rho$pp$beta0
        rho$dpars <- seq_along(rho$pp$theta)
        if (nAGQ > 1L) {
            if (length(vals$reTrms$flist) != 1L || length(vals$reTrms$cnms[[1]]) != 1L)
                stop("nAGQ > 1 is only available for models with a single, scalar random-effects term")
            rho$fac <- vals$reTrms$flist[[1]]
        }
        devfun <- mkdevfun(rho, nAGQ)
        if (devFunOnly) return(devfun)

        opt <- optwrap(optimizer[[2]], devfun, par=c(rho$pp$theta, rho$beta0),
                       lower=rho$lower, control=control,
                       adj=TRUE, verbose=verbose)


    }
    mkMerMod(environment(devfun), opt, vals$reTrms, vals$frame, mc)
}## {nlmer}

## global variables defs to make codetools/R CMD check happier.
## FIXME: does putting globalVariables() stuff here interfere with Roxygen?
## ?globalVariables says that fields and methods in reference classes are
## "handled automatically by ‘setRefClass()’ and friends, using the
##  supplied field and method names" -- perhaps there's a better way to do this?
if (getRversion()<="2.15.0")  {
    ## dummy
    globalVariables <- function(...) {}
}
if(FALSE)## not ok for roxygen2
globalVariables(c("pp","resp","lp0","pwrssUpdate","compDev",
                  "baseOffset","GQmat","fac","nlmerAGQ","tolPwrss",
                  "dpars","verbose"),
                package="lme4")

##' Create a deviance evaluation function from a predictor and a response module
##'
##' From an merMod object create an R function that takes a single argument,
##' which is the new parameter value, and returns the deviance.
##'
##' The function returned by \code{mkdevfun} evaluates the deviance of the model
##' represented by the predictor module, \code{pp}, and the response module,
##' \code{resp}.
##'
##' For \code{\link{lmer}} model objects the argument of the resulting function
##' is the variance component parameter, \code{theta}, with lower bound.  For
##' \code{glmer} or \code{nlmer} model objects with \code{nAGQ = 0} the argument
##' is also \code{theta}.  However, when nAGQ > 0 the argument is \code{c(theta,
##' beta)}.
##'
##' @param rho an environment containing \code{pp}, a prediction module,
##'     typically of class \code{\linkS4class{merPredD}} and \code{resp}, a response
##'     module, e.g., of class \code{\linkS4class{lmerResp}}.
##' @param nAGQ scalar integer - the number of adaptive Gauss-Hermite quadrature
##'     points.  A value of 0 indicates that both the fixed-effects parameters
##'     and the random effects are optimized by the iteratively reweighted least
##'     squares algorithm.
##' @param verbose Logical: print verbose output?
##' @return A function of one numeric argument.
##' @seealso \code{\link{lmer}}, \code{\link{glmer}} and \code{\link{nlmer}}
##' @keywords models
##' @examples
##'
##' (dd <- lmer(Yield ~ 1|Batch, Dyestuff, devFunOnly=TRUE))
##' dd(0.8)
##' minqa::bobyqa(1, dd, 0)
mkdevfun <- function(rho, nAGQ=1L, verbose=0) {
    ## FIXME: should nAGQ be automatically embedded in rho?
    stopifnot(is.environment(rho), is(rho$resp, "lmResp"))

    ## R CMD check  "no visible binding for global variable" ...
    ## MM *preferred* to globalVariables()
    fac <- pp <- resp <- lp0 <- compDev <- dpars <- baseOffset <- tolPwrss <-
	pwrssUpdate <- ## <-- even though it's a function below
	GQmat <- nlmerAGQ <- NULL

    ## The deviance function (to be returned):
    ff <-
    if (is(rho$resp, "lmerResp")) {
	rho$lmer_Deviance <- lmer_Deviance
	function(theta) .Call(lmer_Deviance, pp$ptr(), resp$ptr(), as.double(theta))
    } else if (is(rho$resp, "glmResp")) {
	if (nAGQ == 0L)
	    function(theta) {
		resp$updateMu(lp0)
		pp$setTheta(theta)
		pwrssUpdate(pp, resp, 1e-7, GHrule(0L), compDev, verbose)
	    }
	else
	    function(pars) {
		resp$updateMu(lp0)
		pp$setTheta(as.double(pars[dpars])) # theta is first part of pars
		resp$setOffset(baseOffset + pp$X %*% as.numeric(pars[-dpars]))
		pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose)
	    }
    } else if (is(rho$resp, "nlsResp")) {
	if (nAGQ < 2L) {
	    rho$nlmerLaplace <- nlmerLaplace
	    switch(nAGQ + 1L,
			 function(theta)
			 .Call(nlmerLaplace, pp$ptr(), resp$ptr(), as.double(theta),
			       as.double(u0), beta0, verbose, FALSE, tolPwrss),
			 function(pars)
			 .Call(nlmerLaplace, pp$ptr(), resp$ptr(), pars[dpars], u0,
			       pars[-dpars], verbose, TRUE, tolPwrss))
	} else {
	    rho$nlmerAGQ <- nlmerAGQ
	    rho$GQmat	 <- GHrule(nAGQ)
	    function(pars)
		.Call(nlmerAGQ, pp$ptr(), resp$ptr(), fac, GQmat, pars[dpars],
		      u0, pars[-dpars], tolPwrss)
	}
    }
    else stop("code not yet written")
    environment(ff) <- rho
    ff
}

## Determine a step factor that will reduce the pwrss
##
## The penalized, weighted residual sum of squares (pwrss) is the sum
## of the weighted residual sum of squares from the resp module and
## the squared length of u from the predictor module.  The predictor module
## contains a base value and an increment for the coefficients.
## @title Determine a step factor
## @param pp predictor module
## @param resp response module
## @param verbose logical value determining verbose output
## @return NULL if successful
## @note Typically all this is done in the C++ code.
##     The R code is for debugging and comparisons of
##     results.
## stepFac <- function(pp, resp, verbose, maxSteps = 10) {
##     stopifnot(is.numeric(maxSteps), maxSteps >= 2)
##     pwrss0 <- resp$wrss() + pp$sqrL(0)
##     for (fac in 2^(-(0:maxSteps))) {
## 	wrss <- resp$updateMu(pp$linPred(fac))
## 	pwrss1 <- wrss + pp$sqrL(fac)
## 	if (verbose > 3L)
## 	    cat(sprintf("pwrss0=%10g, diff=%10g, fac=%6.4f\n",
## 			pwrss0, pwrss0 - pwrss1, fac))
## 	if (pwrss1 <= pwrss0) {
## 	    pp$installPars(fac)
## 	    return(NULL)
## 	}
##     }
##     stop("step factor reduced below ",signif(2^(-maxSteps),2)," without reducing pwrss")
## }

RglmerWrkIter <- function(pp, resp, uOnly=FALSE) {
    pp$updateXwts(resp$sqrtWrkWt())
    pp$updateDecomp()
    pp$updateRes(resp$wtWrkResp())
    if (uOnly) pp$solveU() else pp$solve()
    resp$updateMu(pp$linPred(1))	# full increment
    resp$resDev() + pp$sqrL(1)
}

glmerPwrssUpdate <- function(pp, resp, tol, GQmat, compDev=TRUE, grpFac=NULL, verbose=0) {
    nAGQ <- nrow(GQmat)
    if (compDev) {
        if (nAGQ < 2L)
            return(.Call(glmerLaplace, pp$ptr(), resp$ptr(), nAGQ, tol, verbose))
        return(.Call(glmerAGQ, pp$ptr(), resp$ptr(), tol, GQmat, grpFac, verbose))
    }
    oldpdev <- .Machine$double.xmax
    uOnly   <- nAGQ == 0L
    i <- 0
    repeat {
        ## oldu <- pp$delu
        ## olddelb <- pp$delb
        pdev <- RglmerWrkIter(pp, resp, uOnly=uOnly)
        if (verbose>2) cat(i,": ",pdev,"\n",sep="")
        ## check convergence first so small increases don't trigger errors
        if (abs((oldpdev - pdev) / pdev) < tol)
            break
        ## if (pdev > oldpdev) {
        ##     ## try step-halving
        ##     ## browser()
        ##     k <- 0
        ##     while (k < 10 && pdev > oldpdev) {
        ##         pp$setDelu((oldu + pp$delu)/2.)
        ##         if (!uOnly) pp$setDelb((olddelb + pp$delb)/2.)
        ##         pdev <- RglmerWrkIter(pp, resp, uOnly=uOnly)
        ##         k <- k+1
        ##     }
        ## }
        if (pdev>oldpdev) stop("PIRLS update failed")
        oldpdev <- pdev
        i <- i+1
    }
    resp$Laplace(pp$ldL2(), 0., pp$sqrL(1))  ## FIXME: should 0. be pp$ldRX2 ?
}

## create a deviance evaluation function that uses the sigma parameters
## df2 <- function(dd) {
##     stopifnot(is.function(dd),
## 	      length(formals(dd)) == 1L,
## 	      is((rem <- (rho <- environment(dd))$rem), "Rcpp_reModule"),
## 	      is((fem <- rho$fem), "Rcpp_deFeMod"),
## 	      is((resp <- rho$resp), "Rcpp_lmerResp"),
## 	      all((lower <- rem$lower) == 0))
##     Lind <- rem$Lind
##     n <- length(resp$y)
##     function(pars) {
## 	sigma <- pars[1]
## 	sigsq <- sigma * sigma
## 	sigmas <- pars[-1]
## 	theta <- sigmas/sigma
## 	rem$theta <- theta
## 	resp$updateMu(numeric(n))
## 	solveBetaU(rem, fem, resp$sqrtXwt, resp$wtres)
## 	resp$updateMu(rem$linPred1(1) + fem$linPred1(1))
## 	n * log(2*pi*sigsq) + (resp$wrss + rem$sqrLenU)/sigsq + rem$ldL2
##     }
## }

## bootMer() ---> now in ./bootMer.R


## Methods for the merMod class

## Anova for merMod objects
##
## @title anova() for merMod objects
## @param a merMod object
## @param ...	further such objects
## @return an "anova" data frame; the traditional (S3) result of anova()
anovaLmer <- function(object, ...) {
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    modp <- as.logical(.sapply(dots, is, "merMod")) | as.logical(.sapply(dots, is, "lm"))
    if (any(modp)) {			# multiple models - form table
	opts <- dots[!modp]
	mods <- c(list(object), dots[modp])
	## model names
	mNms <- .sapply(as.list(mCall)[c(FALSE, TRUE, modp)], deparse)
	names(mods) <- sub("@env$", '', mNms) # <- hack
	mods <- lapply(mods, refitML)

	devs <- sapply(mods, deviance)
	llks <- lapply(mods, logLik)
	ii <- order(Df <- .sapply(llks, attr, "df"))
	mods <- mods[ii]
	llks <- llks[ii]
	Df   <- Df  [ii]
	calls <- lapply(mods, getCall)
	data <- lapply(calls, "[[", "data")
	if (any(data != data[[1]]))
	    stop("all models must be fit to the same data object")
	header <- paste("Data:", data[[1]])
	subset <- lapply(calls, "[[", "subset")
	if (any(subset != subset[[1]]))
	    stop("all models must use the same subset")
	if (!is.null(subset[[1]]))
	    header <-
		c(header, paste("Subset", deparse(subset[[1]]),
				sep = ": "))
	llk <- unlist(llks)
	chisq <- 2 * pmax(0, c(NA, diff(llk)))
	dfChisq <- c(NA, diff(Df))
	val <- data.frame(Df = Df,
			  AIC = .sapply(llks, AIC),
			  BIC = .sapply(llks, BIC),
                          logLik = llk,
			  deviance = -2*llk,
			  Chisq = chisq,
			  "Chi Df" = dfChisq,
			  "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower.tail = FALSE),
			  row.names = names(mods), check.names = FALSE)
	class(val) <- c("anova", class(val))
	attr(val, "heading") <-
	    c(header, "Models:",
	      paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
				     "[[", "formula"), deparse), length))),
		    unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
		    sep = ": "))
	return(val)
    }
    else { ## ------ single model ---------------------
	dc <- getME(object, "devcomp")
	p <- dc$dims["p"]
        X <- getME(object, "X")
	asgn <- attr(X, "assign")
	stopifnot(length(asgn) == (p <- dc$dims["p"]))
	ss <- as.vector(object@pp$RX() %*% object@beta)^2
	names(ss) <- colnames(X)
	terms <- terms(object)
        nmeffects <- attr(terms, "term.labels")
	if ("(Intercept)" %in% names(ss))
	    nmeffects <- c("(Intercept)", nmeffects)
	ss <- unlist(lapply(split(ss, asgn), sum))
	stopifnot(length(ss) == length(nmeffects))
	df <- unlist(lapply(split(asgn,	 asgn), length))
	## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	ms <- ss/df
	f <- ms/(sigma(object)^2)
	## P <- pf(f, df, dfr, lower.tail = FALSE)
	## table <- data.frame(df, ss, ms, dfr, f, P)
	table <- data.frame(df, ss, ms, f)
	dimnames(table) <-
	    list(nmeffects,
		 ## c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		 c("Df", "Sum Sq", "Mean Sq", "F value"))
	if ("(Intercept)" %in% nmeffects)
	    table <- table[-match("(Intercept)", nmeffects), ]
	attr(table, "heading") <- "Analysis of Variance Table"
	class(table) <- c("anova", "data.frame")
	table
    }
}## {anovaLmer}

##' @importFrom stats anova
##' @S3method anova merMod
anova.merMod <- anovaLmer

##' @S3method as.function merMod
as.function.merMod <- function(x, ...) {
    rho <- list2env(list(resp=x@resp$copy(),
                           pp=x@pp$copy(),
                           beta0=x@beta,
                           u0=x@u), parent=as.environment("package:lme4"))
    mkdevfun(rho, getME(x, "devcomp")$dims["nAGQ"])
}

## coef() method for all kinds of "mer", "*merMod", ... objects
## ------  should work with fixef() + ranef()  alone
coefMer <- function(object, ...)
{
    if (length(list(...)))
	warning('arguments named "', paste(names(list(...)), collapse = ", "),
                '" ignored')
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    val <- lapply(ref, function(x)
		  fef[rep.int(1L, nrow(x)),,drop = FALSE])
    for (i in seq(a = val)) {
	refi <- ref[[i]]
	row.names(val[[i]]) <- row.names(refi)
	nmsi <- colnames(refi)
	if (!all(nmsi %in% names(fef)))
	    stop("unable to align random and fixed effects")
	for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
    }
    class(val) <- "coef.mer"
    val
} ##  {coefMer}

##' @importFrom stats coef
##' @S3method coef merMod
coef.merMod <- coefMer

## FIXME: should these values (i.e. ML criterion for REML models
##  and vice versa) be computed and stored in the object in the first place?
##' @importFrom stats deviance
##' @S3method deviance merMod
deviance.merMod <- function(object, REML = NULL, ...) {
    ## cf. (1) lmerResp::Laplace in respModule.cpp
    ##     (2) section 5.6 of lMMwR, listing lines 34-42
    if (isTRUE(REML) && !isLMM(object))
        stop("can't compute REML deviance for a non-LMM")
    cmp <- object@devcomp$cmp
    if (is.null(REML) || is.na(REML[1]))
        REML <- isREML(object)
    if (REML) {
        if (isREML(object)) {
            cmp["REML"]
        } else {
            lnum <- log(2*pi*(cmp["pwrss"]))
            n <- object@devcomp$dims["n"]
            nmp <- n-length(object@beta)
            unname(cmp["ldL2"]+cmp["ldRX2"]+nmp*(1.+lnum-log(nmp)))
        }
    } else {
        if (!isREML(object)) {
            cmp[["dev"]]
        } else {
            n <- object@devcomp$dims["n"]
            lnum <- log(2*pi*(cmp["pwrss"]))
            unname(cmp["ldL2"]+n*(1+lnum-log(n)))
        }
    }
}

##' @importFrom stats drop1
##' @S3method drop1 merMod
drop1.merMod <- function(object, scope, scale = 0, test = c("none", "Chisq"),
                         k = 2, trace = FALSE, ...) {
    ## FIXME: incorporate na.predict() stuff?
    tl <- attr(terms(object), "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
	if(!is.character(scope)) {
	    scope <- attr(terms(getFixedFormula(update.formula(object, scope))),
                                "term.labels")
        }
	if(!all(match(scope, tl, 0L) > 0L))
	    stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1L, ncol = 2L,
                  dimnames =  list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- nobs(object, use.fallback = TRUE)
    env <- environment(formula(object))
    for(i in seq_along(scope)) {  ## was seq(ns), failed on empty scope
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep='')
	    utils::flush.console()
        }
        nfit <- update(object, as.formula(paste("~ . -", tt)),
                       evaluate = FALSE)
	nfit <- eval(nfit, envir = env) # was  eval.parent(nfit)
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        nnew <- nobs(nfit, use.fallback = TRUE)
        if(all(is.finite(c(n0, nnew))) && nnew != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[1L , 1L] - ans[, 1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        ## BMB: hack to extract safe_pchisq
        P[nas] <- stats:::safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    } else if (test == "F") {
        ## FIXME: allow this if denominator df are specified externally?
        stop("F test STUB -- unfinished maybe forever")
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        ## BMB: hack to extract safe_pchisq
        P[nas] <- stats:::safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(F)")] <- list(dev, P)
    }
    head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

##' @importFrom stats extractAIC
##' @S3method extractAIC merMod
extractAIC.merMod <- function(fit, scale = 0, k = 2, ...) {
    L <- logLik(refitML(fit))
    edf <- attr(L,"df")
    c(edf,-2*L + k*edf)
}

##' @importFrom stats family
##' @S3method family merMod
family.merMod <- function(object, ...) family(object@resp, ...)

##' @S3method family glmResp
family.glmResp <- function(object, ...) object$family

##' @S3method family lmResp
family.lmResp <- function(object, ...) gaussian()

##' @S3method family nlsResp
family.nlsResp <- function(object, ...) gaussian()

##' @importFrom stats fitted
##' @S3method fitted merMod
fitted.merMod <- function(object, ...) object@resp$mu

##' Extract the fixed-effects estimates
##'
##' Extract the estimates of the fixed-effects parameters from a fitted model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixed.effects fixef.merMod
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a named, numeric vector of fixed-effects estimates.
##' @keywords models
##' @examples
##' fixef(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' @importFrom nlme fixef
##' @export fixef
##' @method fixef merMod
##' @export
fixef.merMod <- function(object, ...)
    structure(object@beta, names = dimnames(object@pp$X)[[2]])

getFixedFormula <- function(form) {
    form[[3]] <- if (is.null(nb <- nobars(form[[3]]))) 1 else nb
    form
}

##' @importFrom stats formula
##' @S3method formula merMod
formula.merMod <- function(x, fixed.only=FALSE, ...) {
    form <- formula(getCall(x),...)
    if (fixed.only) {
        form <- getFixedFormula(form)
    }
    form
}

##' @S3method isREML merMod
isREML.merMod <- function(x, ...) as.logical(x@devcomp$dims["REML"])

##' @S3method isGLMM merMod
isGLMM.merMod <- function(x,...) {
  as.logical(x@devcomp$dims["GLMM"])
  ## or: is(x@resp,"glmResp")
}

##' @S3method isNLMM merMod
isNLMM.merMod <- function(x,...) {
  as.logical(x@devcomp$dims["NLMM"])
  ## or: is(x@resp,"nlsResp")
}

##' @S3method isLMM merMod
isLMM.merMod <- function(x,...) {
  !isGLMM(x) && !isNLMM(x)
  ## or: is(x@resp,"lmerResp") ?
}

##' @importFrom stats logLik
##' @S3method logLik merMod
logLik.merMod <- function(object, REML = NULL, ...) {
    if (is.null(REML) || is.na(REML[1]))
        REML <- isREML(object)
    val <- -deviance(object, REML = REML)/2
    dc <- object@devcomp
    dims <- dc$dims
    attr(val, "nall") <- attr(val, "nobs") <- nrow(object@frame) ## FIXME use nobs() ?
    attr(val, "df") <- length(object@beta) + length(object@theta) + dims[["useSc"]]
    class(val) <- "logLik"
    val
}

##' @importFrom stats logLik
##' @S3method model.frame merMod
model.frame.merMod <- function(formula, ...) formula@frame

##' @importFrom stats model.matrix
##' @S3method model.matrix merMod
model.matrix.merMod <- function(object, ...) object@pp$X

##' @importFrom stats nobs
##' @S3method nobs merMod
nobs.merMod <- function(object, ...) nrow(object@frame)

##' @importFrom nlme ranef
##' @export ranef
NULL

##' Extract the modes of the random effects
##'
##' A generic function to extract the conditional modes of the random effects
##' from a fitted model object.  For linear mixed models the conditional modes
##' of the random effects are also the conditional means.
##'
##' If grouping factor i has k levels and j random effects per level the ith
##' component of the list returned by \code{ranef} is a data frame with k rows
##' and j columns.  If \code{postVar} is \code{TRUE} the \code{"postVar"}
##' attribute is an array of dimension j by j by k.  The kth face of this array
##' is a positive definite symmetric j by j matrix.  If there is only one
##' grouping factor in the model the variance-covariance matrix for the entire
##' random effects vector, conditional on the estimates of the model parameters
##' and on the data will be block diagonal and this j by j matrix is the kth
##' diagonal block.  With multiple grouping factors the faces of the
##' \code{"postVar"} attributes are still the diagonal blocks of this
##' conditional variance-covariance matrix but the matrix itself is no longer
##' block diagonal.
##' @name ranef
##' @aliases ranef ranef.merMod
##' @param object an object of a class of fitted models with random effects,
##' typically an \code{"\linkS4class{merMod}"} object.
##' @param postVar an optional logical argument indicating if the conditional
##' variance-covariance matrices, also called the \dQuote{posterior variances},
##' of the random effects should be added as an attribute.  Default is
##' \code{FALSE}.
##' @param drop an optional logical argument indicating components of the return
##' value that would be data frames with a single column, usually a column
##' called \sQuote{\code{(Intercept)}}, should be returned as named vectors.
##' @param whichel an optional character vector of names of grouping factors for
##' which the random effects should be returned.  Defaults to all the grouping
##' factors.
##' @param \dots some methods for this generic function require additional
##' arguments.
##' @return A list of data frames, one for each grouping factor for the random
##' effects.  The number of rows in the data frame is the number of levels of
##' the grouping factor.  The number of columns is the dimension of the random
##' effect associated with each level of the factor.
##'
##' If \code{postVar} is \code{TRUE} each of the data frames has an attribute
##' called \code{"postVar"} which is a three-dimensional array with symmetric
##' faces.
##'
##' When \code{drop} is \code{TRUE} any components that would be data frames of
##' a single column are converted to named numeric vectors.
##' @note To produce a \dQuote{caterpillar plot} of the random effects apply
##' \code{\link[lattice:xyplot]{dotplot}} to the result of a call to
##' \code{ranef} with \code{postVar = TRUE}.
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
##' fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
##' ranef(fm1)
##' str(rr1 <- ranef(fm1, postVar = TRUE))
##' dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
##' if(FALSE) { ##-- postVar=TRUE is not yet implemented for multiple terms -- FIXME
##' str(ranef(fm2, postVar = TRUE))
##' }
##' op <- options(digits = 4)
##' ranef(fm3, drop = TRUE)
##' options(op)
##' @keywords models methods
##' @method ranef merMod
##' @export
ranef.merMod <- function(object, postVar = FALSE, drop = FALSE,
			 whichel = names(ans), ...)
{
    ans <- object@pp$b(1.)
    if (!is.null(object@flist)) {
	## evaluate the list of matrices
	levs <- lapply(fl <- object@flist, levels)
	asgn <- attr(fl, "assign")
	cnms <- object@cnms
	nc <- sapply(cnms, length)
	nb <- nc * (nl <- unlist(lapply(levs, length))[asgn])
	nbseq <- rep.int(seq_along(nb), nb)
	ml <- split(ans, nbseq)
	for (i in seq_along(ml))
	    ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE,
			      dimnames = list(NULL, cnms[[i]]))
	## create a list of data frames corresponding to factors
	ans <- lapply(seq_along(fl),
		      function(i)
		      data.frame(do.call(cbind, ml[asgn == i]),
				 row.names = levs[[i]],
				 check.names = FALSE))
	names(ans) <- names(fl)
					# process whichel
	stopifnot(is(whichel, "character"))
	whchL <- names(ans) %in% whichel
	ans <- ans[whchL]

	if (postVar) {
            sigsqr <- sigma(object)^2
	    vv <- .Call(merPredDcondVar, object@pp$ptr(), as.environment(rePos$new(object)))
	    for (i in names(ans)) ## seq_along(ans))
                attr(ans[[i]], "postVar") <- vv[[i]] * sigsqr
	}
	if (drop)
	    ans <- lapply(ans, function(el)
		      {
			  if (ncol(el) > 1) return(el)
			  pv <- drop(attr(el, "postVar"))
			  el <- drop(as.matrix(el))
			  if (!is.null(pv))
			      attr(el, "postVar") <- pv
			  el
		      })
	class(ans) <- "ranef.mer"
    }
    ans
}## ranef.merMod


##' @method refit merMod
##' @rdname refit
##' @export
refit.merMod <- function(object, newresp=NULL, ...)
{
    rr <- object@resp$copy()

    if (!is.null(newresp)) {

        if (!is.null(na.act <- attr(object@frame,"na.action"))) {
            ## will only get here if na.action is 'na.omit' or 'na.exclude'
            if (is.matrix(newresp)) {
                newresp <- newresp[-na.act,]
            } else newresp <- newresp[-na.act]
        }

        if (isGLMM(object) && rr$family$family=="binomial") {
            if (is.matrix(newresp) && ncol(newresp)==2) {
                ntot <- rowSums(newresp)
                ## FIXME: test what happens for (0,0) columns
                newresp <- newresp[,1]/ntot
                rr$setWeights(ntot)
            }
            if (is.factor(newresp)) {
                ## FIXME: would be better to do this consistently with
                ## whatever machinery is used in glm/glm.fit/glmer ... ??
                newresp <- as.numeric(newresp)-1
            }
        }
        stopifnot(length(newresp <- as.numeric(as.vector(newresp))) == length(rr$y))
        rr$setResp(newresp)
    }

    pp        <- object@pp$copy()
    dc        <- object@devcomp
    nAGQ      <- dc$dims["nAGQ"]
    nth       <- dc$dims["nth"]
    verbose   <- list(...)$verbose
    if (is.null(verbose)) verbose <- 0L
    devlist <- list(pp=pp, resp=rr, u0=pp$u0, verbose=verbose, dpars=seq_len(nth))
    if (isGLMM(object)) {
        baseOffset <- object@resp$offset
        devlist <- c(list(tolPwrss=unname(dc$cmp["tolPwrss"]),
                          compDev=unname(dc$dims["compDev"]),
                          nAGQ=unname(nAGQ),
                          lp0=object@resp$eta - baseOffset,
                          baseOffset=baseOffset,
                          pwrssUpdate=glmerPwrssUpdate,
                          ## save GQmat in the object and use that instead of nAGQ
                          GQmat=GHrule(nAGQ)), devlist)
    }
    ff <- mkdevfun(list2env(devlist),nAGQ=nAGQ)
    xst       <- rep.int(0.1, nth)
    x0        <- pp$theta
    lower     <- object@lower
    if (!is.na(nAGQ) && nAGQ > 0L) {
        xst   <- c(xst, sqrt(diag(pp$unsc())))
        x0    <- c(x0, unname(fixef(object)))
        lower <- c(lower, rep(-Inf,length(x0)-length(lower)))
    }
    control <- list(...)$control
    if (is.null(control)) control <- list()
    control <- c(control,list(xst=0.2*xst, xt=xst*0.0001))
    ## FIXME: generic optimizer stuff
### FIXME: Probably should save the control settings and the optimizer name in the merMod object
    opt <- Nelder_Mead(ff, x0, lower=lower, control=control)
    if (isGLMM(object)) rr$setOffset(baseOffset)
    mkMerMod(environment(ff), opt,
             list(flist=object@flist, cnms=object@cnms, Gp=object@Gp, lower=object@lower),
             object@frame, getCall(object))
}

##-- BUG in roxygen2: If we use  @S3method instead of @method,
##-- the \usage{ ... } will have
##-- refitML.merMod(..) instead of \method{refitML}{mermod}(..)
##' @param optimizer a string indicating the optimizer to be used.
##' @method refitML merMod
##' @rdname refitML
##' @export
refitML.merMod <- function (x, optimizer="bobyqa", ...) {
    ## FIXME: optimizer is set to 'bobyqa' for back-compatibility, but that's not
    ##  consistent with lmer (default NM).  Should be based on internally stored 'optimizer' value
    if (!isREML(x)) return(x)
    stopifnot(is(rr <- x@resp, "lmerResp"))
    rho <- new.env(parent=parent.env(environment()))
    rho$resp <- new(class(rr), y=rr$y, offset=rr$offset, weights=rr$weights, REML=0L)
    xpp <- x@pp
    rho$pp <- new(class(xpp), X=xpp$X, Zt=xpp$Zt, Lambdat=xpp$Lambdat,
                  Lind=xpp$Lind, theta=xpp$theta, n=nrow(xpp$X))
    devfun <- mkdevfun(rho, 0L)
    opt <- optwrap(optimizer, devfun, x@theta, lower=x@lower)
    ##  opt <- bobyqa(x@theta, devfun, x@lower)
    n <- length(rr$y)
    pp <- rho$pp
    p <- ncol(pp$X)
    dims <- c(N=n, n=n, nmp=n-p, nth=length(pp$theta), p=p, q=nrow(pp$Zt),
	      nAGQ=NA_integer_, useSc=1L, reTrms=length(x@cnms),
	      spFe=0L, REML=0L, GLMM=0L, NLMM=0L)
    wrss <- rho$resp$wrss()
    ussq <- pp$sqrL(1)
    pwrss <- wrss + ussq
    cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss, ussq=ussq,
	     pwrss=pwrss, drsum=NA, dev=opt$fval, REML=NA,
	     sigmaML=sqrt(pwrss/n), sigmaREML=sqrt(pwrss/(n-p)))
### FIXME: Should modify the call slot to set REML=FALSE.  It is
### tricky to do so without causing the call to be evaluated
    new("lmerMod", call=x@call, frame=x@frame, flist=x@flist,
	cnms=x@cnms, theta=pp$theta, beta=pp$delb, u=pp$delu,
	lower=x@lower, devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=rho$resp)
}

##' @importFrom stats residuals
##' @S3method residuals merMod
residuals.merMod <-
    function(object,
             type = c("deviance", "pearson", "working", "response", "partial"),
             ...)
        residuals(object@resp, match.arg(type), ...)

##' @S3method residuals lmResp
residuals.lmResp <- function(object, type = c("deviance", "pearson",
                                     "working", "response", "partial"),
                             ...) {
### FIXME: This should be extended with na.resid but need to store na.action
    switch(match.arg(type),
           working =,
           response = object$y - object$mu,
           deviance =,
           pearson = object$wtres,
           partial = .NotYetImplemented())
}

## FIXME: document somewhere that residuals(glmerfit) returns deviance residuals by default
##' @S3method residuals glmResp
residuals.glmResp <- function(object, type = c("deviance", "pearson",
                                      "working", "response", "partial"),
                              ...) {
    type <- match.arg(type)
    y <- object$y
    r <- object$wtres
    mu <- object$mu
    switch(type,
           deviance = {
               d.res <- sqrt(object$devResid())
               ifelse(y > mu, d.res, -d.res)
           },
           pearson = object$wtres,
           working = object$wrkResids(),
           response = y - mu,
           partial = .NotYetImplemented())
}

##' @S3method sigma merMod
sigma.merMod <- function(object, ...) {
    dc <- object@devcomp
    dd <- dc$dims
    if(dd[["useSc"]])
        dc$cmp[[ifelse(dd[["REML"]], "sigmaREML", "sigmaML")]] else 1.
}

##' @importFrom stats simulate
NULL
##' Simulate responses from the model represented by a fitted model object
##'
##' @title Simulate responses from a \code{\linkS4class{merMod}} object
##' @param object a fitted model object
##' @param nsim positive integer scalar - the number of responses to simulate
##' @param seed an optional seed to be used in \code{set.seed} immediately
##'     before the simulation so as to generate a reproducible sample.
##' @param use.u (logical) generate new random-effects values (FALSE) or
##'     generate a simulation condition on the current random-effects estimates (TRUE)?
##' @param ... optional additional arguments, none are used at present
##' @examples
##' ## test whether fitted models are consistent with the
##' ##  observed number of zeros in CBPP data set:
##' gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'              data = cbpp, family = binomial)
##' gg <- simulate(gm1,1000)
##' zeros <- sapply(gg,function(x) sum(x[,"incidence"]==0))
##' plot(table(zeros))
##' abline(v=sum(cbpp$incidence==0),col=2)
##' @method simulate merMod
##' @export
simulate.merMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE, ...) {
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
	      is(object, "merMod"))
	      ## i.e. not yet for glmer etc:
    ## is(object@resp, "lmerResp"))
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
	runif(1) # initialize the RNG if necessary
    RNGstate <- .Random.seed

    sigma <- sigma(object)
    n <- nrow(X <- getME(object, "X"))
    if (is.null(nm <- names(fitted(object)))) nm <- seq(n)
    # fixed-effect contribution
    etasim.fix <- as.vector(X %*% getME(object, "beta"))
    if (length(offset <- getME(object,"offset"))>0) {
      etasim.fix <- etasim.fix+offset
    }
    etasim.reff <- ## UNSCALED random-effects contribution:
      if(use.u) {
        getME(object, "u")
      } else {
        U <- getME(object, "Z") %*% getME(object, "Lambda")
        q <- ncol(U)
        as(U %*% matrix(rnorm(q * nsim), ncol = nsim), "matrix")
      }
    if (is(object@resp,"lmerResp")) {
      ## result will be matrix  n x nsim :
      val <- etasim.fix + sigma * (etasim.reff +
        ## residual contribution:
        matrix(rnorm(n * nsim), ncol = nsim))
    } else if (is(object@resp,"glmResp")) {
      ## GLMM
      ## n.b. DON'T scale random-effects (???)
      	      etasim <- etasim.fix+etasim.reff
	      family <- object@call$family
	      if(is.symbol(family)) family <- as.character(family)
	      if(is.character(family))
		  family <- get(family, mode = "function", envir = parent.frame(2))
	      if(is.function(family)) family <- family()
              if(is.language(family)) family <- eval(family)
	      if(is.null(family$family)) stop("'family' not recognized")
	      musim <- family$linkinv(etasim)
	      ntot <- length(musim) ## FIXME: or could be dims["n"]?
              ## FIXME: is it possible to leverage family$simulate ... ???
              val <- switch(family$family,
			    poisson=rpois(ntot,lambda=musim),
			    binomial={
                              w <- weights(object)
                              Y <- rbinom(ntot,prob=musim,size=w)
                              resp <- model.response(object@frame)
                              if (!is.matrix(resp)) {  ## bernoulli, or weights specified
                                if (is.factor(resp)) {
                                  if (any(weights(object)!=1)) stop("non-uniform weights with factor response??")
                                  factor(levels(resp)[Y+1],levels=levels(resp))
                                } else {
                                  Y/w
                                }
                              } else {
                                ## FIXME: should "N-size" (column 2) be named?
                                ## copying structures from stats/R/family.R
                                nresp <- nrow(resp)
                                YY <- cbind(Y, w - Y)
                                yy <- lapply(split(YY,gl(nsim,nresp,2*nsim*nresp)),
                                             matrix, ncol=2,
                                             dimnames=list(NULL,colnames(resp)))
                                names(yy) <- paste("sim",seq_along(yy),sep="_")
                                yy
                              }
                            },
			    stop("simulation not implemented for family",
				 family$family))
            } else {
              stop("simulate method for NLMMs not yet implemented")
            }
    ## from src/library/stats/R/lm.R
    if(!is.list(val)) {
      dim(val) <- c(n, nsim)
      val <- as.data.frame(val)
    }
    else
      class(val) <- "data.frame"
    names(val) <- paste("sim", seq_len(nsim), sep="_")
    row.names(val) <- nm
    attr(val, "seed") <- RNGstate
    val
  }

##' @importFrom stats terms
##' @S3method terms merMod
terms.merMod <- function(x, fixed.only=TRUE, ...) {
  if (fixed.only) {
      terms.formula(formula(x,fixed.only=TRUE))
  } else attr(x@frame,"terms")
}

##' @importFrom stats update
##' @S3method update merMod
update.merMod <- function(object, formula., ..., evaluate = TRUE) {
    if (is.null(call <- getCall(object)))
	stop("object should contain a 'call' component")
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

## This is modeled a bit after	print.summary.lm :
## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
			correlation = NULL, symbolic.cor = FALSE,
			signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    cat(sprintf("%s ['%s']\n",so$methTitle, class(x)))
    if (!is.null(f <- so$family)) {
	cat(" Family:", f)
        if (!(is.null(ll <- so$link))) cat(" (", ll, ")")
        cat("\n")
    }
    if (!is.null(cc <- so$call$formula))
	cat("Formula:", deparse(cc),"\n")
    ## if (!is.null(so$family)) {
    ##     cat("Family: ",so$family,
    ##         " (link=",so$link,")\n",
    ##         sep="")
    ## }
    if (!is.null(cc <- so$call$data))
	cat("   Data:", deparse(cc), "\n")
    if (!is.null(cc <- so$call$subset))
	cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
    cat("\n")
    tab <- so$AICtab
    if (length(tab) == 1 && names(tab) == "REML")
	cat("REML criterion at convergence:", round(tab, 4), "\n")
    else print(round(so$AICtab, 4))
    cat("\nRandom effects:\n")
    print(formatVC(so$varcor, digits = digits, useScale = so$useScale),
	  quote = FALSE, digits = digits, ...)

    ngrps <- so$ngrps
    cat(sprintf("Number of obs: %d, groups: ", so$devcomp$dims[["n"]]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    p <- nrow(so$coefficients)
    if (p > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so$coefficients, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(!is.logical(correlation)) { # default
	    correlation <- p <= 20
	    if(!correlation) {
		nam <- deparse(substitute(x)) # << TODO: improve if this is called from show()
		cat(sprintf(paste("\nCorrelation matrix not shown by default, as p = %d > 20.",
				  "Use print(%s, correlation=TRUE)  or",
				  "    vcov(%s)	 if you need it\n", sep="\n"),
			    p, nam, nam))
	    }
	}
	if(correlation) {
	    if(is.null(VC <- so$vcov)) VC <- vcov(x)
	    corF <- VC@factors$correlation
	    if (is.null(corF)) {
		cat("\nCorrelation of Fixed Effets is not available\n")
	    }
	    else {
		p <- ncol(corF)
		if (p > 1) {
		    rn <- rownames(so$coefficients)
		    rns <- abbreviate(rn, minlength=11)
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			corf <- as(corF, "matrix")
			dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			print(symnum(corf))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       ncol = p,
				       dimnames = list(rns, abbreviate(rn, minlength=6)))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}## printMerenv()

##' @S3method print merMod
print.merMod <- printMerenv

##' @exportMethod show
setMethod("show",  "merMod", function(object) printMerenv(object))

##' @S3method print summary.mer
print.summary.mer <- printMerenv

##' Return the deviance component list
##'
##' A fitted model of class \code{\linkS4class{merMod}} has a \code{devcomp}
##' slot as described in the value section.
##' @title Extract the deviance component list
##' @param x a fitted model of class \code{\linkS4class{merMod}}
##' @return a list with components
##' \item{dims}{a named integer vector of various dimensions}
##' \item{cmp}{a named numeric vector of components of the deviance}
##' @export
##' @note This function is deprecated, use \code{getME(., "devcomp")}
devcomp <- function(x) {
    .Deprecated("getME(., \"devcomp\")")
    stopifnot(is(x, "merMod"))
    x@devcomp
}

##' @exportMethod getL
setMethod("getL", "merMod", function(x) {
    .Deprecated("getME(., \"L\")")
    getME(x, "L")
})

##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##'
##' Extract (or \dQuote{get}) \dQuote{components} -- in a generalized sense --
##' from a fitted mixed-effects model, i.e. (in this version of the package)
##' from an object of class \code{"\linkS4class{merMod}"}.
##'
##' The goal is to provide \dQuote{everything a user may want} from a fitted
##' \code{"merMod"} object \emph{as far} as it is not available by methods, such
##' as \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{vcov}}, etc.
##'
##' @aliases getME getL getL,merMod-method
##' @param object a fitted mixed-effects model of class
##' \code{"\linkS4class{merMod}"}, i.e. typically the result of
##' \code{\link{lmer}()}, \code{\link{glmer}()} or \code{\link{nlmer}()}.
##' @param name a character string specifying the name of the
##' \dQuote{component}.  Possible values are:\cr
##' \describe{
##'     \item{X}{fixed-effects model matrix}
##'     \item{Z}{random-effects model matrix}
##'     \item{Zt}{transpose of random-effects model matrix}
##'     \item{u}{conditional mode of the \dQuote{spherical} random effects variable}
##'     \item{Gp}{groups pointer vector.  A pointer to the beginning of each group
##'               of random effects corresponding to the random-effects terms.}
##'     \item{L}{sparse Cholesky factor of the penalized random-effects model.}
##'     \item{Lambda}{relative covariance factor of the random effects.}
##'     \item{Lambdat}{transpose of the relative covariance factor of the random effects.}
##'     \item{Lind}{index vector for inserting elements of \eqn{\theta}{theta} into the
##'                 nonzeros of \eqn{\Lambda}{Lambda}}
##'     \item{A}{Scaled sparse model matrix (class
##'      \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
##'      the unit, orthogonal random effects, \eqn{U},
##'       equal to \code{getME(.,"Zt") \%*\% getME(.,"Lambdat")}}
##'     \item{RX}{Cholesky factor for the fixed-effects parameters}
##'     \item{RZX}{cross-term in the full Cholesky factor}
##'     \item{flist}{a list of the grouping variables (factors) involved in the random effect terms}
##'     \item{beta}{fixed-effects parameter estimates (identical to the result of \code{\link{fixef}}, but without names)}
##'     \item{theta}{random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term}
##'     \item{n_rtrms}{number of random-effects terms}
##'     \item{is_REML}{same as the result of \code{\link{isREML}}}
##'     \item{devcomp}{a list consisting of a named numeric vector, \dQuote{cmp}, and
##'                    a named integer vector, \dQuote{dims}, describing the fitted model}
##'     \item{offset}{model offset}
##'     \item{lower}{lower bounds on model parameters (random effects parameters only)}
##' }
##' @return Unspecified, as very much depending on the \code{\link{name}}.
##' @seealso \code{\link{getCall}()},
##' More standard methods for mer objects, such as \code{\link{ranef}},
##' \code{\link{fixef}}, \code{\link{vcov}}, etc.:
##' see \code{methods(class="merMod")}
##' @keywords utilities
##' @examples
##'
##' ## shows many methods you should consider *before* using getME():
##' methods(class = "merMod")
##'
##' (fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
##' Z <- getME(fm1, "Z")
##' stopifnot(is(Z, "CsparseMatrix"),
##'           c(180,36) == dim(Z),
##' 	  all.equal(fixef(fm1), getME(fm1, "beta"),
##' 		    check.attr=FALSE, tol = 0))
##'
##' ## All that can be accessed [potentially ..]:
##' (nmME <- eval(formals(getME)$name))
##' \dontshow{
##' ## internal consistency check ensuring that all work:
##' ## "try(.)" because some are not yet implemented:
##' str(parts <- sapply(nmME, function(nm) try(getME(fm1, nm)),
##'                     simplify=FALSE))
##' }% dont..
##'
##' @export
getME <- function(object,
		  name = c("X", "Z","Zt", "u",
		  "Gp",
		  "L", "Lambda", "Lambdat", "Lind", "A",
		  "RX", "RZX",
                  "flist",
                  "beta", "theta",
		  "REML", "n_rtrms", "is_REML", "devcomp",
                    "offset", "lower"))
{
    if(missing(name)) stop("'name' must not be missing")
    stopifnot(length(name <- as.character(name)) == 1,
	      is(object, "merMod"))
    name <- match.arg(name)
    rsp  <- object@resp
    PR   <- object@pp
    dc   <- object@devcomp
    cmp  <- dc $ cmp
    dims <- dc $ dims
    switch(name,
	   "X" = PR$X, ## ok ? - check -- use model.matrix() method instead?
	   "Z" = t(PR$Zt),
	   "Zt"= PR$Zt,
           "u" = object@u,
	   "L"= PR$ L(),
	   "Lambda"= t(PR$ Lambdat),
	   "Lambdat"= PR$ Lambdat,
           "A" = PR$Lambdat %*% PR$Zt,
           "Lind" = PR$ Lind,
	   "RX" = PR $ RX(), ## FIXME - add the column names and row names, either in the C++ or the R method
	   "RZX" = PR $ RZX, ## FIXME - add column names

           "Gp" = object@Gp,
           "flist" = object@flist,
	   "beta" = object@beta,
           "theta"= {
               tt <- object@theta
               nc <- c(unlist(mapply(function(g,e) {
                   mm <- outer(e,e,paste,sep=".")
                   diag(mm) <- e
                   mm <- mm[lower.tri(mm,diag=TRUE)]
                   paste(g,mm,sep=".")
               },
                                     names(object@cnms),object@cnms)))
               names(tt) <- nc
               tt
           }, ## *OR*  PR $ theta  --- which one ??  Should be the same.

	   "REML" = dims["REML"],
	   "is_REML" = isREML(object),

	   "n_rtrms" = length(object@flist),  ## should this be length(object@cnms) instead?

           ## Yes, assuming that you want the number of random-effects
           ## terms in the formula.  Multiple terms with the same
           ## grouping factors are allowed.

           "devcomp" = dc,
           "offset" = rsp$offset,
           "lower" = object@lower,
            ## FIXME: current version gives lower bounds for theta parameters only -- these must be extended for [GN]LMMs -- give extended value including -Inf values for beta values?
	   "..foo.." =# placeholder!
	   stop(gettextf("'%s' is not implemented yet",
			 sprintf("getME(*, \"%s\")", name))),
	   ## otherwise
	   stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
			name, class(object))))
}## {getME}

##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod
##' @importClassesFrom Matrix dgCMatrix dpoMatrix corMatrix
NULL

## Extract the conditional variance-covariance matrix of the fixed-effects
## parameters
##
## @title Extract conditional covariance matrix of fixed effects
## @param sigma numeric scalar, the residual standard error
## @param unsc matrix of class \code{"\linkS4class{dpoMatrix}"}, the
##     unscaled variance-covariance matrix
## @param nmsX character vector of column names of the model matrix
## @param correlation logical scalar, should the correlation matrix
##     also be evaluated.
## @param ... additional, optional parameters.  None are used at present.
mkVcov <- function(sigma, unsc, nmsX, correlation = TRUE, ...) {
    V <- sigma^2 * unsc
    if(is.null(rr <- tryCatch(as(V, "dpoMatrix"),
			      error = function(e) NULL)))
	stop("Computed variance-covariance matrix is not positive definite")
    dimnames(rr) <- list(nmsX, nmsX)
    if(correlation)
	rr@factors$correlation <- as(rr, "corMatrix")
    rr
}

##' @importFrom stats vcov
##' @S3method vcov merMod
vcov.merMod <- function(object, correlation = TRUE, sigm = sigma(object), ...)
    mkVcov(sigm, unsc = object@pp$unsc(), nmsX = colnames(object@pp$X),
	   correlation=correlation, ...)

##' @importFrom stats vcov
##' @S3method vcov summary.mer
vcov.summary.mer <- function(object, correlation = TRUE, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of merMod object")
    object$vcov
}

mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    ans <- lapply(ncseq, function(i)
	      {
		  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		  Li <- diag(nrow = nc[i])
		  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		  rownames(Li) <- cnms[[i]]
		  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
		  val <- tcrossprod(sc * Li) # variance-covariance
		  stddev <- sqrt(diag(val))
		  correl <- t(val / stddev)/stddev
		  diag(correl) <- 1
		  attr(val, "stddev") <- stddev
		  attr(val, "correlation") <- correl
		  val
	      })
    if(is.character(nms)) names(ans) <- nms
    attr(ans, "sc") <- sc
    ans
}

##' Extract variance and correlation components
##'
##' This function calculates the estimated variances, standard deviations, and
##' correlations between the random-effects terms in a mixed-effects model, of
##' class \code{\linkS4class{merMod}} (linear, generalized or nonlinear).  The
##' within-group error variance and standard deviation are also calculated.
##'
##' @name VarCorr
##' @aliases VarCorr VarCorr.merMod
##' @param x a fitted model object, usually an object inheriting from class
##' \code{\linkS4class{merMod}}.
##' @param sigma an optional numeric value used as a multiplier for the standard
##' deviations.  Default is \code{1}.
##' @param rdig an optional integer value specifying the number of digits used
##' to represent correlation estimates.  Default is \code{3}.
##' @return a list of matrices, one for each random effects grouping term.
##' For each grouping term, the standard deviations and correlation matrices for each grouping term
##' are stored as attributes \code{"stddev"} and \code{"correlation"}, respectively, of the
##' variance-covariance matrix, and
##' the residual standard deviation is stored as attribute \code{"sc"}
##' (for \code{glmer} fits, this attribute stores the scale parameter of the model).
##' @author This is modeled after \code{\link[nlme]{VarCorr}} from package
##' \pkg{nlme}, by Jose Pinheiro and Douglas Bates.
##' @seealso \code{\link{lmer}}, \code{\link{nlmer}}
##' @examples
##' data(Orthodont, package="nlme")
##' fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
##' VarCorr(fm1)
##' @keywords models
##' @importFrom nlme VarCorr
##' @export VarCorr
##' @method VarCorr merMod
##' @export
VarCorr.merMod <- function(x, sigma, rdig)# <- 3 args from nlme
{
  ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    if (is.null(cnms <- x@cnms))
	stop("VarCorr methods require reTrms, not just reModule")
    if(missing(sigma)) # "bug": fails via default 'sigma=sigma(x)'
	sigma <- lme4::sigma(x)  ## FIXME: do we still need lme4:: ?
    nc <- sapply(cnms, length)	  # no. of columns per term
    m <- mkVarCorr(sigma, cnms=cnms, nc=nc, theta = x@theta,
	      nms = {fl <- x@flist; names(fl)[attr(fl, "assign")]})
    attr(m,"useSc") <- as.logical(x@devcomp$dims["useSc"])
    class(m) <- "VarCorr.merMod"
    m
}

## FIXME: should ... go to formatVC or to print ... ?
##' @S3method print VarCorr.merMod
print.VarCorr.merMod <- function(x,digits = max(3, getOption("digits") - 2), ...) {
  print(formatVC(x, digits = digits, useScale = attr(x,"useSc"),  ...),quote=FALSE)
}

## Compute standard errors of fixed effects from an merMod object
##
## @title Standard errors of fixed effects
## @param object "merMod" object,
## @param ... additional, optional arguments.  None are used at present.
## @return numeric vector of length length(fixef(.))
unscaledVar <- function(object, ...) {
    stopifnot(is(object, "merMod"))
    sigma(object) * diag(object@pp$unsc())
}

### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
		     useScale) {
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), if(useScale) list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }))[, -maxlen, drop = FALSE]
	if (nrow(corr) < nrow(reMat))
	    corr <- rbind(corr, matrix("", nrow = nrow(reMat) - nrow(corr), ncol = ncol(corr)))
	colnames(corr) <- rep.int("", ncol(corr))
	colnames(corr)[1] <- "Corr"
	cbind(reMat, corr)
    } else reMat
}

##' @S3method summary merMod
summary.merMod <- function(object, ...)
{
    resp <- object@resp
    devC <- object@devcomp
    dd <- devC$dims
    cmp <- devC$cmp
    useSc <- as.logical(dd["useSc"])
    sig <- sigma(object)
    REML <- isREML(object)

    link <- fam <- NULL
    if(is(resp, "glmResp")) {
        fam <- resp$family$family
        link <- resp$family$link
    }
    coefs <- cbind("Estimate" = fixef(object),
		   "Std. Error" = sig * sqrt(diag(object@pp$unsc())))
    if (nrow(coefs) > 0) {
	coefs <- cbind(coefs, coefs[,1]/coefs[,2], deparse.level=0)
	colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
        if (isGLMM(object)) {
          coefs <- cbind(coefs,2*pnorm(abs(coefs[,3]),lower.tail=FALSE))
          colnames(coefs)[4] <- c("Pr(>|z|)")
        }
    }
    mName <- paste(switch(1L + dd["GLMM"] * 2L + dd["NLMM"],
			  "Linear", "Nonlinear",
			  "Generalized linear", "Generalized nonlinear"),
		   "mixed model fit by",
		   ifelse(REML, "REML", "maximum likelihood"))
    llik <- logLik(object)   # returns NA for a REML fit - maybe change?
    AICstats <- {
	if (REML) cmp["REML"] # do *not* show likelihood stats here
	else {
	    c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
	      deviance = deviance(object))
	}
    }
    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
					# use S3 class for now
    structure(list(methTitle=mName, devcomp=devC,
                   isLmer=is(resp, "lmerResp"), useScale=useSc,
		   logLik=llik, family=fam, link=link,
		   ngrps=sapply(object@flist, function(x) length(levels(x))),
		   coefficients=coefs, sigma=sig,
		   vcov=vcov(object, correlation=TRUE, sigm=sig),
		   varcor=varcor, # and use formatVC(.) for printing.
		   AICtab=AICstats, call=object@call
		   ), class = "summary.mer")
}

##' @S3method summary summary.mer
summary.summary.mer <- function(object, varcov = FALSE, ...) {
    if(varcov && is.null(object$vcov))
	object$vcov <- vcov(object, correlation=TRUE, sigm = object$sigma)
    object
}

### Plots for the ranef.mer class ----------------------------------------

##' @importFrom lattice dotplot
##' @S3method  dotplot ranef.mer
dotplot.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
	if (is.null(se)) return(list())
	x <- as.numeric(x)
	hw <- 1.96 * as.numeric(se[subscripts])
	list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16,
			 horizontal = TRUE, col = dot.symbol$col,
			 lty = dot.line$lty, lwd = dot.line$lwd,
			 col.line = dot.line$col, levels.fos = unique(y),
			 groups = NULL, ...)
    {
	x <- as.numeric(x)
	y <- as.numeric(y)
	dot.line <- trellis.par.get("dot.line")
	dot.symbol <- trellis.par.get("dot.symbol")
	sup.symbol <- trellis.par.get("superpose.symbol")
	panel.abline(h = levels.fos, col = col.line, lty = lty, lwd = lwd)
	panel.abline(v = 0, col = col.line, lty = lty, lwd = lwd)
	if (!is.null(se)) {
	    se <- as.numeric(se[subscripts])
	    panel.segments( x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
	}
	panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x, ...) {
	ss <- stack(x)
	ss$ind <- factor(as.character(ss$ind), levels = colnames(x))
	ss$.nn <- rep.int(reorder(factor(rownames(x)), x[[1]]), ncol(x))
	se <- NULL
	if (!is.null(pv <- attr(x, "postVar")))
	    se <- unlist(lapply(1:(dim(pv)[1]), function(i) sqrt(pv[i, i, ])))
	dotplot(.nn ~ values | ind, ss, se = se,
		prepanel = prepanel.ci, panel = panel.ci,
		xlab = NULL, ...)
    }
    lapply(x, f, ...)
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
qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
	x <- as.numeric(x)
	se <- as.numeric(se[subscripts])
	hw <- 1.96 * se
	list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
	panel.grid(h = -1,v = -1)
	panel.abline(v = 0)
	x <- as.numeric(x)
	y <- as.numeric(y)
	se <- as.numeric(se[subscripts])
	panel.segments(x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
	panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x) {
	if (!is.null(pv <- attr(x, "postVar"))) {
	    cols <- 1:(dim(pv)[1])
	    se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
	    nr <- nrow(x)
	    nc <- ncol(x)
	    ord <- unlist(lapply(x, order)) +
		rep((0:(nc - 1)) * nr, each = nr)
	    rr <- 1:nr
	    ind <- gl(ncol(x), nrow(x), labels = names(x))
	    xyplot(rep(qnorm((rr - 0.5)/nr), ncol(x)) ~ unlist(x)[ord] | ind[ord],
		   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
		   scales = list(x = list(relation = "free")),
		   ylab = "Standard normal quantiles",
		   xlab = NULL, ...)
	} else {
	    qqmath(~values|ind, stack(x),
		   scales = list(y = list(relation = "free")),
		   xlab = "Standard normal quantiles",
		   ylab = NULL, ...)
	}
    }
    lapply(x, f)
}

##' @importFrom graphics plot
##' @S3method plot coef.mer
plot.coef.mer <- function(x, y, ...)
{
    ## remove non-varying columns from frames
    reduced <- lapply(x, function(el)
		      el[, !sapply(el, function(cc) all(cc == cc[1]))])
    plot.ranef.mer(reduced, ...)
}

##' @importFrom lattice dotplot
##' @S3method dotplot coef.mer
dotplot.coef.mer <- function(x, data, ...) {
    mc <- match.call()
    mc[[1]] <- as.name("dotplot.ranef.mer")
    eval(mc)
}

##' @importFrom stats weights
##' @S3method weights merMod
weights.merMod <- function(object, ...) {
  object@resp$weights
}


getOptfun <- function(optimizer) {
  if (is.character(optimizer)) {
    optfun <- try(get(optimizer),silent=TRUE)
  } else optfun <- optimizer
  if (inherits(optfun,"try-error")) stop("couldn't find optimizer function ",optimizer)
  if (!is(optfun,"function")) stop("non-function specified as optimizer")
  if (any(is.na(match(c("fn","par","lower","control"),names(formals(optfun))))))
    stop("optimizer function must use (at least) formal parameters 'fn', 'par', 'lower', 'control'")
  optfun
}

optwrap <- function(optimizer, fn, par, lower=-Inf, upper=Inf,
                    control=list(), adj=FALSE, verbose=0L) {
    ## control must be specified if adj==TRUE;
    ##  otherwise this is a fairly simple wrapper
    optfun <- getOptfun(optimizer)

    lower <- rep(lower, length.out=length(par))
    upper <- rep(upper, length.out=length(par))

    if (adj && is.character(optimizer))
        ## control parameter tweaks: only for second round in nlmer, glmer
        switch(optimizer,
               bobyqa = {
                   if(!is.numeric(control$rhobeg)) control$rhobeg <- 0.0002
                   if(!is.numeric(control$rhoend)) control$rhoend <- 2e-7
               },
               Nelder_Mead = {
                   if (is.null(control$xst))
                       xst <- c(rep.int(0.1, length(environment(fn)$pp$theta)),  ## theta parameters
                                sqrt(diag(environment(fn)$pp$unsc())))
                   control$xst <- 0.2*xst
                   control$verbose <- verbose
                   if (is.null(control$xt)) control$xt <- control$xst*5e-4
               })
    if (optimizer=="bobyqa" && all(par==0)) par[] <- 0.001  ## minor kluge
    arglist <- list(fn=fn, par=par, lower=lower, upper=upper, control=control)
    ## optimx: must pass method in control (?) because 'method' was previously
    ## used in lme4 to specify REML vs ML
    ## FIXME: test -- does deparse(substitute(...)) clause work?
    if (optimizer=="optimx" || deparse(substitute(optimizer))=="optimx") {
        if (is.null(method <- control$method))
            stop("must specify 'method' explicitly for optimx")
        arglist$control$method <- NULL
        arglist <- c(arglist,list(method=method))
    }
    ## FIXME: test!  effects of multiple warnings??
    ## may not need to catch warnings after all??
    curWarning <- NULL
    opt <- withCallingHandlers(do.call(optfun,arglist),
                               warning = function(w) {
                                   ## browser()
                                   cat("caught warning:",w$message,"\n")
                                   assign("curWarning",w$message,parent.frame())
                                   invokeRestart("muffleWarning")
                               })
    ## if (!is.null(curWarning)) browser()
    ## FIXME: set code to warn on convergence !=0
    ## post-fit tweaking
    if (optimizer=="bobyqa") {
        opt$convergence <- opt$ierr
    }
    if (optimizer=="optimx") {
        optr <- lapply(opt,"[[",1)[c("par","fvalues","conv")]
        optr$message <- attr(opt,"details")[[1]]$message
        opt <- optr
    }
    if (opt$conv!=0) {
        wmsg <- paste("convergence code",opt$conv,"from",optimizer)
        if (!is.null(opt$msg)) wmsg <- paste0(wmsg,": ",opt$msg)
        warning(wmsg)
    }
    attr(opt,"control") <- control
    opt
}


