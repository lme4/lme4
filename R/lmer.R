## FIXME: documentation still refers to \linkS4class quite a bit, inappropriately
## FIXME: need to document S3 methods better (can we pull from r-forge version?)

##'
##' Fit a linear mixed model (LMM)
##'
##' @title Fit Linear Mixed-Effects Models
##' @concept LMM
##' @aliases lmer
##' @param formula a two-sided linear formula object describing both the fixed-effects and
##'    fixed-effects part of the model, with the response on the left of a
##'    \code{~} operator and the terms, separated by \code{+} operators, on
##'    the right.  Random-effects terms are distinguished by vertical bars
##'    (\code{"|"}) separating expressions for design matrices from
##'    grouping factors.
##' @param data an optional data frame containing the variables named in
##'    \code{formula}.  By default the variables are taken from the environment
##'    from which \code{lmer} is called. While \code{data} is optional,
##'    the package authors \emph{strongly} recommend its use,
##'    especially when later applying methods such as
##'    \code{update} and \code{drop1} to the fitted model
##' (\emph{such methods are not guaranteed to work properly if \code{data} is omitted}).
##' If \code{data} is omitted, variables will be taken from the environment
##' of \code{formula} (if specified as a formula) or from the parent frame
##' (if specified as a character vector).
##' @param REML logical scalar - Should the estimates be chosen to optimize
##'    the REML criterion (as opposed to the log-likelihood)?
##' @param control a list (of correct class, resulting from
##'    \code{\link{lmerControl}()} or \code{\link{glmerControl}()} respectively)
##'    containing control parameters, including the nonlinear optimizer to be
##'    used and parameters to be passed through to the nonlinear optimizer, see
##'    the \code{*lmerControl} documentation for details.
##' @param start a named list of starting values for the parameters in the
##'    model.  For \code{lmer} this can be a numeric vector or a list with one
##'    component named \code{"theta"}.
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
##'     data contain \code{NA}s.  The default action (\code{na.omit},
##'     inherited from the 'factory fresh' value of \code{getOption("na.action")})
##'     strips any observations with any missing values in any variables.
##' @param offset this can be used to specify an \emph{a priori} known component
##'     to be included in the linear predictor during fitting. This should be
##'     \code{NULL} or a numeric vector of length equal to the number of cases.
##'     One or more \code{\link{offset}} terms can be included in the formula
##'     instead or as well, and if more than one is specified their sum is used.
##'     See \code{\link{model.offset}}.
##' @param contrasts an optional list. See the \code{contrasts.arg} of
##'     \code{model.matrix.default}.
##' @param devFunOnly logical - return only the deviance evaluation function. Note that
##'    because the deviance function operates on variables stored in its environment,
##'    it may not return \emph{exactly} the same values on subsequent calls (but the results should always be within machine tolerance).
##' @param \dots other potential arguments.  A \code{method} argument was used
##'    in earlier versions of the package. Its functionality has been replaced by
##'    the \code{REML} argument.
##' @return An object of class \code{merMod}, for which many
##'    methods are available (e.g. \code{methods(class="merMod")})
##' @seealso \code{\link[stats]{lm}}
##' @keywords models
##' @details
##' \itemize{
##' \item{If the \code{formula} argument is specified as a character vector,
##' the function will attempt to coerce it to a formula. However, this is
##' not recommended (users who want to construct formulas by pasting together
##' components are advised to use \code{\link{as.formula}} or \code{\link{reformulate}}); model fits will
##' work but subsequent methods such as \code{\link{drop1}}, \code{\link{update}}
##' may fail.}
##' \item{Unlike some simpler modeling frameworks such as \code{\link{lm}}
##' and \code{\link{glm}} which automatically detect perfectly collinear
##' predictor variables, \code{[gn]lmer} cannot handle design matrices of
##' less than full rank.  For example, in cases of models with interactions
##' that have unobserved combinations of levels, it is up to the user to
##' define a new variable (for example creating
##' \code{ab} within the data from the results of \code{interaction(a,b,drop=TRUE)}).
##' }
##' \item{the deviance function returned when \code{devFunOnly} is \code{TRUE}
##' takes a single numeric vector argument, representing the \code{theta} vector.
##' This vector defines the scaled variance-covariance matrices of the random effects,
##' in the Cholesky parameterization.  For models with only simple (intercept-only) random effects,
##' \code{theta} is a vector of the standard deviations of the random effects. For more
##' complex or multiple random effects, running \code{getME(.,"theta")} to
##' retrieve the \code{theta} vector for a fitted model and examining the
##' names of the vector is probably the easiest way to determine the correspondence
##' between the elements of the \code{theta} vector and elements of the lower
##' triangles of the Cholesky factors of the random effects.}
##' }
##' @examples
##' ## linear mixed models - reference values from older code
##' (fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
##' fm1_ML <- update(fm1,REML=FALSE)
##' (fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' anova(fm1, fm2)
##' @export
##' @importFrom minqa bobyqa
lmer <- function(formula, data=NULL, REML = TRUE,
                 control = lmerControl(), start = NULL,
                 verbose = 0L, subset, weights, na.action, offset,
                 contrasts = NULL, devFunOnly=FALSE,
                 ...)
{

    mc <- mcout <- match.call()
    missCtrl <- missing(control)
    ## see functions in modular.R for the body ...
    if (!missCtrl && !inherits(control, "lmerControl")) {
        if(!is.list(control)) stop("'control' is not a list; use lmerControl()")
        ## back-compatibility kluge
	warning("passing control as list is deprecated: please use lmerControl() instead",
		immediate.=TRUE)
        control <- do.call(lmerControl, control)
    }
    if (!is.null(list(...)[["family"]])) {
       warning("calling lmer with 'family' is deprecated; please use glmer() instead")
       mc[[1]] <- quote(lme4::glmer)
       if(missCtrl) mc$control <- glmerControl()
       return(eval(mc, parent.frame(1L)))
    }
    mc$control <- control ## update for  back-compatibility kluge

    ## https://github.com/lme4/lme4/issues/50
    ## parse data and formula
    mc[[1]] <- quote(lme4::lFormula)
    lmod <- eval(mc, parent.frame(1L))  ## parse data and formula
    mcout$formula <- lmod$formula
    lmod$formula <- NULL

    ## create deviance function for covariance parameters (theta)
    devfun <- do.call(mkLmerDevfun,
		      c(lmod,
			list(start=start,verbose=verbose,control=control)))
    if (devFunOnly) return(devfun)
    ## optimize deviance function over covariance parameters
    opt <- optimizeLmer(devfun,
                        optimizer=control$optimizer,
                        restart_edge=control$restart_edge,
                        control=control$optCtrl,
                        verbose=verbose,
                        start=start)
    mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, mcout) ## prepare output
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
##'
##' @title Fit Generalized Linear Mixed-Effects Models
##' @concept GLMM
##' @param family a GLM family, see \code{\link[stats]{glm}} and
##'    \code{\link[stats]{family}}.
##' @param nAGQ integer scalar - the number of points per axis for evaluating
##'    the adaptive Gauss-Hermite approximation to the log-likelihood.
##'    Defaults to 1, corresponding to the Laplace
##'    approximation.  Values greater than 1 produce greater accuracy in
##'    the evaluation of the log-likelihood at the expense of speed.  A value
##'    of zero uses a faster but less exact form of parameter estimation for
##'    GLMMs by optimizing the random effects and the fixed-effects coefficients
##'    in the penalized iteratively reweighted least squares step.
##' @param start a named list of starting values for the parameters in the
##'    model, or a numeric vector. A numeric \code{start} argument
##'    will be used as the starting value of \code{theta}.  If \code{start}
##'    is a list, the \code{theta} element (a numeric vector) is used
##'    as the starting value for the first optimization step (default=1
##'    for diagonal elements and 0 for off-diagonal elements of the
##'    lower Cholesky factor); the fitted value of \code{theta} from
##'    the first step, plus \code{start[["fixef"]]},
##'    are used as starting values for the second optimization step.
##'    If \code{start} has both \code{fixef} and \code{theta}
##'    elements, the first optimization step is skipped. For more details
##'    or finer control of optimization, see \code{\link{modular}}.
##' @param mustart optional starting values on the scale of the conditional mean,
##'    as in \code{\link[stats]{glm}}; see there for details.
##' @param etastart optional starting values on the scale of the unbounded
##'    predictor as in \code{\link[stats]{glm}}; see there for details.
##' @param \dots other potential arguments.  A \code{method} argument was used
##'    in earlier versions of the package. Its functionality has been replaced by
##'    the \code{nAGQ} argument.
##' @inheritParams lmer
##' @return An object of class \code{glmerMod}, for which many
##'    methods are available (e.g. \code{methods(class="glmerMod")})
##' @seealso \code{\link{lmer}} (for details on formulas and parameterization); \code{\link[stats]{glm}}
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
##'
##' ## glmer and glm log-likelihoods are consistent
##' gm1Devfun <- update(gm1,devFunOnly=TRUE)
##' gm0 <- glm(cbind(incidence, size - incidence) ~ period,
##'            family = binomial, data = cbpp)
##' ## evaluate GLMM deviance at RE variance=theta=0, beta=(GLM coeffs)
##' gm1Dev0 <- gm1Devfun(c(0,coef(gm0)))
##' ## compare
##' stopifnot(all.equal(gm1Dev0,c(-2*logLik(gm0))))
##'
##' @export
glmer <- function(formula, data=NULL, family = gaussian,
                  control = glmerControl(), start = NULL, verbose = 0L, nAGQ = 1L,
                  subset, weights, na.action, offset,
                  contrasts = NULL, mustart, etastart, devFunOnly = FALSE, ...)
{
    if (!inherits(control, "glmerControl")) {
	if(!is.list(control)) stop("'control' is not a list; use glmerControl()")
	## back-compatibility kluge
	msg <- "Use control=glmerControl(..) instead of passing a list"
	if(length(cl <- class(control))) msg <- paste(msg, "of class", dQuote(cl[1]))
	warning(msg, immediate.=TRUE)
	control <- do.call(glmerControl, control)
    }
    mc <- mcout <- match.call()

    ## family-checking code duplicated here and in glFormula (for now) since
    ## we really need to redirect at this point; eventually deprecate formally
    ## and clean up
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if( is.function(family)) family <- family()
    if (isTRUE(all.equal(family, gaussian()))) {
        ## redirect to lmer (with warning)
        warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
                " please call lmer() directly")
        mc[[1]] <- quote(lme4::lmer)
        mc["family"] <- NULL            # to avoid an infinite loop
        return(eval(mc, parent.frame()))
    }

    ## see https://github.com/lme4/lme4/issues/50
    ## parse the formula and data
    mc[[1]] <- quote(lme4::glFormula)
    glmod <- eval(mc, parent.frame(1L))
    mcout$formula <- glmod$formula
    glmod$formula <- NULL

    ## create deviance function for covariance parameters (theta)

    devfun <- do.call(mkGlmerDevfun, c(glmod, list(verbose=verbose,
                                                   control=control,
                                                   nAGQ = 0)))
    if (nAGQ==0 && devFunOnly) return(devfun)
    ## optimize deviance function over covariance parameters

    if (is.list(start) && !is.null(start$fixef))
        if (nAGQ==0) stop("should not specify both start$fixef and nAGQ==0")

    opt <- optimizeGlmer(devfun,
                         optimizer = control$optimizer[[1]],
                         restart_edge=control$restart_edge,
                         control = control$optCtrl,
                         start=start,
                         nAGQ = 0,
                         verbose=verbose)

    if(nAGQ > 0L) {

        start <- updateStart(start,theta=opt$par)

        # update deviance function to include fixed effects as inputs
        devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)

        if (devFunOnly) return(devfun)
        # reoptimize deviance function over covariance parameters and fixed effects
        opt <- optimizeGlmer(devfun,
                             optimizer = control$optimizer[[2]],
                             restart_edge=control$restart_edge,
                             control = control$optCtrl,
                             start=start,
                             nAGQ=nAGQ,
                             verbose = verbose,
                             stage=2)
    }
    # prepare output
    mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr, mcout)

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
##' @note Adaptive Gauss-Hermite quadrature (\code{nAGQ}>1) is not currently implemented for \code{nlmer}.
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
nlmer <- function(formula, data=NULL, control = nlmerControl(), start = NULL, verbose = 0L,
                  nAGQ = 1L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly = FALSE, ...)
{

    vals <- nlformula(mc <- match.call())
    p <- ncol(X <- vals$X)
    if ((rankX <- rankMatrix(X)) < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))

    rho <- list2env(list(verbose=verbose,
                         tolPwrss=0.001, # this is reset to the tolPwrss argument's value later
                         resp=vals$resp,
                         lower=vals$reTrms$lower),
                    parent=parent.frame())
    rho$pp <- do.call(merPredD$new,
                      c(vals$reTrms[c("Zt","theta","Lambdat","Lind")],
                        list(X=X, n=length(vals$respMod$mu), Xwts=vals$respMod$sqrtXwt,
                             beta0=qr.coef(qr(X), unlist(lapply(vals$pnames, get,
                             envir = rho$resp$nlenv))))))
    rho$u0 <- rho$pp$u0
    rho$beta0 <- rho$pp$beta0
    devfun <- mkdevfun(rho, 0L, verbose, control) # deviance as a function of theta only
    if (devFunOnly && !nAGQ) return(devfun)
    devfun(rho$pp$theta) # initial coarse evaluation to get u0 and beta0
    rho$u0 <- rho$pp$u0
    rho$beta0 <- rho$pp$beta0
    rho$tolPwrss <- control$tolPwrss # Reset control parameter (the initial optimization is coarse)

    opt <- optwrap(control$optimizer[[1]], devfun, rho$pp$theta, rho$lower,
                   control=control$optCtrl, adj=FALSE)
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
        devfun <- mkdevfun(rho, nAGQ, verbose, control)
        if (devFunOnly) return(devfun)

        opt <- optwrap(control$optimizer[[2]], devfun, par=c(rho$pp$theta, rho$beta0),
                       lower=rho$lower, control=control$optCtrl,
                       adj=TRUE, verbose=verbose)


    }
    mkMerMod(environment(devfun), opt, vals$reTrms, vals$frame, mc)
}## {nlmer}

## R 3.1.0 devel [2013-08-05]: This does not help yet
if(getRversion() >= "3.1.0") utils::suppressForeignCheck("nlmerAGQ")
if(getRversion() < "3.1.0") dontCheck <- identity

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
##' @param control list of control parameters, a subset of those specified
##'   by \code{\link{lmerControl}} (\code{tolPwrss} and \code{compDev} for GLMMs,
##' \code{tolPwrss} for NLMMs)
##' @return A function of one numeric argument.
##' @seealso \code{\link{lmer}}, \code{\link{glmer}} and \code{\link{nlmer}}
##' @keywords models
##' @examples
##'
##' (dd <- lmer(Yield ~ 1|Batch, Dyestuff, devFunOnly=TRUE))
##' dd(0.8)
##' minqa::bobyqa(1, dd, 0)
mkdevfun <- function(rho, nAGQ=1L, verbose=0, control=list()) {
    ## FIXME: should nAGQ be automatically embedded in rho?
    stopifnot(is.environment(rho), is(rho$resp, "lmResp"))

    ## silence R CMD check warnings *locally* in this function
    ## (clearly preferred to using globalVariables() !]
    fac <- pp <- resp <- lp0 <- compDev <- dpars <- baseOffset <- tolPwrss <-
	pwrssUpdate <- ## <-- even though it's a function below
	GQmat <- nlmerAGQ <- NULL

    ## The deviance function (to be returned):
    ff <-
    if (is(rho$resp, "lmerResp")) {
	rho$lmer_Deviance <- lmer_Deviance
	function(theta) .Call(lmer_Deviance, pp$ptr(), resp$ptr(), as.double(theta))
    } else if (is(rho$resp, "glmResp")) {
        ## control values will override rho values *if present*
        if (!is.null(tp <- control$tolPwrss)) rho$tolPwrss <- tp
        if (!is.null(cd <- control$compDev)) rho$compDev <- cd
	if (nAGQ == 0L)
	    function(theta) {
		resp$updateMu(lp0)
		pp$setTheta(theta)
		p <- pwrssUpdate(pp, resp, tolPwrss, GHrule(0L),
                            compDev, verbose)
                resp$updateWts()
                p

            }
	else
	    function(pars) {
                ## pp$setDelu(rep(0, length(pp$delu)))
                resp$setOffset(baseOffset)
		resp$updateMu(lp0)
		pp$setTheta(as.double(pars[dpars])) # theta is first part of pars
                spars <- as.numeric(pars[-dpars])
                offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
		resp$setOffset(offset)
		p <- pwrssUpdate(pp, resp, tolPwrss, GQmat,
                            compDev, fac, verbose)
                resp$updateWts()
                p
	    }
    } else if (is(rho$resp, "nlsResp")) {
	if (nAGQ < 2L) {
	    rho$nlmerLaplace <- nlmerLaplace
            rho$tolPwrss <- control$tolPwrss
	    switch(nAGQ + 1L,
			 function(theta)
			 .Call(nlmerLaplace, pp$ptr(), resp$ptr(), as.double(theta),
			       as.double(u0), beta0, verbose, FALSE, tolPwrss),
			 function(pars)
			 .Call(nlmerLaplace, pp$ptr(), resp$ptr(), pars[dpars], u0,
			       pars[-dpars], verbose, TRUE, tolPwrss))
	} else {
            stop("AGQ>1 not yet implemented for nlmer models")
	    rho$nlmerAGQ <- nlmerAGQ
	    rho$GQmat	 <- GHrule(nAGQ)
	    ## function(pars) {
            ## .Call(nlmerAGQ, ## <- dontCheck(nlmerAGQ)  should work according to docs but does not
            ## pp$ptr(), resp$ptr(), fac, GQmat, pars[dpars],
            ## u0, pars[-dpars], tolPwrss)
            ##}
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
        if (is.na(pdev)) stop("encountered NA in PWRSS update")
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
        ## HACK to try to identify model names in situations such as
        ## 'do.call(anova,list(model1,model2))' where the model names
        ## are lost in the call stack ... this doesn't quite work but might
        ## be useful for future attempts?
        ## maxdepth <- -2
        ## depth <- -1
        ## while (depth>=maxdepth &
        ##        all(grepl("S4 object of class structure",mNms))) {
        ##     xCall <- match.call(call=sys.call(depth))
        ##     mNms <- .sapply(as.list(xCall)[c(FALSE, TRUE, modp)], deparse)
        ##     depth <- depth-1
        ## }
        ## if (depth<maxdepth) {
        if (any(duplicated(mNms))) {
            warning("failed to find unique model names, assigning generic names")
            mNms <- paste0("MODEL",seq_along(mNms))
        }
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
	df <- vapply(split(asgn, asgn), length, 1L)
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
    ## FIXME: extract verbose and control
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
    ## check for variables in RE but missing from FE, fill in zeros in FE accordingly
    refnames <- unlist(lapply(ref,colnames))
    nmiss <- length(missnames <- setdiff(refnames,names(fef)))
    if (nmiss >0) {
        fillvars <- setNames(data.frame(rbind(rep(0,nmiss))),missnames)
        fef <- cbind(fillvars,fef)
    }
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
            ## adjust ML results to REML
            lnum <- log(2*pi*(cmp["pwrss"]))
            n <- object@devcomp$dims["n"]
            nmp <- n-length(object@beta)
            unname(cmp["ldL2"]+cmp["ldRX2"]+nmp*(1.+lnum-log(nmp)))
        }
    } else {
        if (!isREML(object)) {
            cmp[["dev"]]
        } else {
            ## adjust REML results to ML
            n <- object@devcomp$dims["n"]
            lnum <- log(2*pi*(cmp["pwrss"]))
            unname(cmp["ldL2"]+n*(1+lnum-log(n)))
        }
    }
}

## copied from stats:::safe_pchisq
safe_pchisq <- function (q, df, ...) {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
}

##' @importFrom stats drop1
##' @S3method drop1 merMod
drop1.merMod <- function(object, scope, scale = 0, test = c("none", "Chisq"),
                         k = 2, trace = FALSE, evalhack="formulaenv", ...) {
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
    env <- environment(formula(object)) # perhaps here is where trouble begins??
    for(i in seq_along(scope)) {  ## was seq(ns), failed on empty scope
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep='')
	    utils::flush.console()
        }
        ## FIXME: make this more robust, somehow?
        ## three choices explored so far:
        ##  (1) evaluate nfit in parent frame: tests in inst/tests/test-formulaEval.R
        ##      will fail on lapply(m_data_List,drop1)
        ##      (formula environment contains r,x,y,z but not d)
        ##  (2) evaluate nfit in frame of formula: tests will fail when data specified and formula is character
        ##  (3) update with data=NULL: fails when ...
        ##
        if (evalhack %in% c("parent","formulaenv")) {
            nfit <- update(object,
                           as.formula(paste("~ . -", tt)),
                           evaluate = FALSE)
            ## nfit <- eval(nfit, envir = env) # was  eval.parent(nfit)
            if (evalhack=="parent") {
                nfit <- eval.parent(nfit)
            } else if (evalhack=="formulaenv") {
                nfit <- eval(nfit,envir=env)
            }
        } else {
            nfit <- update(object,
                           as.formula(paste("~ . -", tt)),data=NULL,
                           evaluate = FALSE)
            nfit <- eval(nfit,envir=env)
        }
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
        ## reconstruct deviance from AIC (ugh)
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        ## BMB: hack to extract safe_pchisq
        P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    } else if (test == "F") {
        ## FIXME: allow this if denominator df are specified externally?
        stop("F test STUB -- unfinished maybe forever")
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        ## BMB: hack to extract safe_pchisq
        P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
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
    if (is.null(form <- attr(x@frame,"formula"))) {
        if (!grepl("lmer$",deparse(getCall(x)[[1]])))
            stop("can't find formula stored in model frame or call")
        form <- as.formula(formula(getCall(x),...))
    }
    if (fixed.only) {
        form <- getFixedFormula(form)
    }
    form
}

##' @S3method isREML merMod
isREML.merMod <- function(x, ...) as.logical(x@devcomp$dims["REML"])

##' @S3method isGLMM merMod
isGLMM.merMod <- function(x,...) {
  as.logical(x@devcomp$dims[["GLMM"]])
  ## or: is(x@resp,"glmResp")
}

##' @S3method isNLMM merMod
isNLMM.merMod <- function(x,...) {
  as.logical(x@devcomp$dims[["NLMM"]])
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
    nobs <- nrow(object@frame) ## FIXME use nobs() ?
    structure(val,
	      nobs = nobs,
	      nall = nobs,
	      df = length(object@beta) + length(object@theta) + dims[["useSc"]],
	      class = "logLik")
}

stripwhite <- function(x) gsub("(^ +| +$)","",x)
##' @importFrom stats logLik
##' @S3method model.frame merMod
model.frame.merMod <- function(formula, fixed.only=FALSE, ...) {
    fr <- formula@frame
    if (fixed.only) {
        ff <- formula(formula,fixed.only=TRUE)
        ## thanks to Thomas Leeper and Roman Luštrik, Stack Overflow
        vars <- rownames(attr(terms.formula(ff), "factors"))
        fr <- fr[vars]
    }
    fr
}

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
##' and j columns.  If \code{condVar} is \code{TRUE} the \code{"postVar"}
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
##' @param condVar an optional logical argument indicating if the conditional
##' variance-covariance matrices of the random effects should be added as an attribute.
##' @param postVar a (deprecated) synonym for \code{condVar}
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
##' If \code{condVar} is \code{TRUE} each of the data frames has an attribute
##' called \code{"postVar"} which is a three-dimensional array with symmetric
##' faces.
##'
##' When \code{drop} is \code{TRUE} any components that would be data frames of
##' a single column are converted to named numeric vectors.
##' @note To produce a \dQuote{caterpillar plot} of the random effects apply
##' \code{\link[lattice:xyplot]{dotplot}} to the result of a call to
##' \code{ranef} with \code{condVar = TRUE}.
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
##' fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
##' ranef(fm1)
##' str(rr1 <- ranef(fm1, condVar = TRUE))
##' dotplot(rr1)  ## default
##' ## specify free scales in order to make Day effects more visible
##' dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
##' if(FALSE) { ##-- condVar=TRUE is not yet implemented for multiple terms -- FIXME
##' str(ranef(fm2, condVar = TRUE))
##' }
##' op <- options(digits = 4)
##' ranef(fm3, drop = TRUE)
##' options(op)
##' @keywords models methods
##' @method ranef merMod
##' @export
ranef.merMod <- function(object, condVar = FALSE, drop = FALSE,
			 whichel = names(ans), postVar = FALSE, ...)
{
    if (!missing(postVar) && missing(condVar)) {
        warning(sQuote("postVar")," is deprecated: please use ",
                sQuote("condVar")," instead")
        condVar <- postVar
    }
    ans <- object@pp$b(1.)
    if (!is.null(object@flist)) {
	## evaluate the list of matrices
	levs <- lapply(fl <- object@flist, levels)
	asgn <- attr(fl, "assign")
	cnms <- object@cnms
	nc <- vapply(cnms, length, 1L)
	nb <- nc * (nl <- vapply(levs, length, 1L)[asgn])
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

	if (condVar) {
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

    ## for backward compatibility/functioning of refit(fit,simulate(fit))
    if (is.list(newresp)) {
        if (length(newresp)==1) {
            newresp <- newresp[[1]]
        } else {
            stop("refit not implemented for lists with length>1: ",
                 "consider ",sQuote("lapply(object,refit)"))
        }
    }
    
    if (!is.null(newresp)) {
        if (!is.null(na.act <- attr(object@frame,"na.action"))) {
            ## will only get here if na.action is 'na.omit' or 'na.exclude'
            if (is.matrix(newresp)) {
                newresp <- newresp[-na.act,]
            } else newresp <- newresp[-na.act]
        }

        if (isGLMM(object) && rr$family$family=="binomial") {
            ## re-do conversion of two-column matrix and factor
            ##  responses to proportion/weights format
            if (is.matrix(newresp) && ncol(newresp)==2) {
                ntot <- rowSums(newresp)
                ## FIXME: test what happens for (0,0) rows
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
    ff <- mkdevfun(list2env(devlist),nAGQ=nAGQ, verbose)
    xst       <- rep.int(0.1, nth)
    x0        <- pp$theta
    lower     <- object@lower
    if (!is.na(nAGQ) && nAGQ > 0L) {
        xst   <- c(xst, sqrt(diag(pp$unsc())))
        x0    <- c(x0, unname(fixef(object)))
        lower <- c(lower, rep(-Inf,length(x0)-length(lower)))
    }
    control <- object@optinfo$control
    newControl <- list(...)$control
    if (!is.null(newControl)) {
        for (i in names(newControl)) {
            control[[i]] <- newControl[[i]]
        }
    }
    ## control <- c(control,list(xst=0.2*xst, xt=xst*0.0001))
    opt <- optwrap(object@optinfo$optimizer,
                   ff, x0, lower=lower, control=control)
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
    xpp <- x@pp$copy()
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

##' residuals of merMod objects
##' @importFrom stats residuals
##' @S3method residuals merMod
##' @method residuals merMod
##' @param object a fitted [g]lmer (\code{merMod}) object
##' @param type type of residuals
##' @param scaled scale residuals by residual standard deviation (=scale parameter)?
##' @param \dots additional arguments (ignored: for method compatibility)
##' @details
##' \itemize{
##' \item The default residual type
##' varies between \code{lmerMod} and \code{glmerMod} objects: they try to
##' mimic \code{\link{residuals.lm}} and \code{\link{residuals.glm}} respectively.
##' In particular, the default \code{type} is \code{"response"},
##' i.e. (observed-fitted) for \code{lmerMod} objects
##' vs. \code{"deviance"} for \code{glmerMod} objects.  \code{type="partial"}
##' is not yet implemented for either type.
##' \item Note that the meaning of \code{"pearson"} residuals differs between
##' \code{\link{residuals.lm}} and \code{\link{residuals.lme}}.  The former
##' returns values scaled by the square root of user-specified weights (if any),
##' but \emph{not} by the residual standard deviation,
##' while the latter returns values scaled by the estimated standard deviation
##' (which will include the effects of any variance structure specified in
##' the \code{weights} argument).  To replicate \code{lme} behaviour, use
##' \code{type="pearson"}, \code{scaled=TRUE}.
##' }
residuals.merMod <-
    function(object,
             type=if (isGLMM(object)) "deviance" else "response",
             scaled=FALSE,
             ...) {
        r <- residuals(object@resp, type,...)
        if (scaled) r <- r/sigma(object)
        if (!is.null(na.action <- attr(model.frame(object),"na.action")))
            r <- naresid(na.action,r)
        r
    }

##' @rdname residuals.merMod
##' @S3method residuals lmResp
##' @method residuals lmResp
residuals.lmResp <- function(object,
                             type = c("working", "response", "deviance",
                             "pearson", "partial"),
                             ...) {
    y <- object$y
    r <- object$wtres
    mu <- object$mu
    switch(match.arg(type),
                    working =,
                    response = y-mu,
                    deviance =,
                    pearson = r,
                    partial = .NotYetImplemented())
}

##' @rdname residuals.merMod
##' @S3method residuals glmResp
##' @method residuals glmResp
residuals.glmResp <- function(object, type = c("deviance", "pearson",
                                      "working", "response", "partial"),
                              ...) {
    type <- match.arg(type)
    y <- object$y
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
        dc$cmp[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]] else 1.
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
##' @param use.u (logical) if \code{TRUE}, generate a simulation conditional on the current
##' random-effects estimates; if \code{FALSE} generate new Normally distributed random-effects values
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
    U <- getME(object, "Z") %*% getME(object, "Lambda")
    u <- if (use.u) {
        rep(getME(object, "u"), nsim)/sigma  ## ??? u is 'spherized' but not scaled ???
    } else {
        rnorm(ncol(U)*nsim)
    }
    etasim.reff <- ## UNSCALED random-effects contribution:
        as(U %*% matrix(u, ncol = nsim), "matrix")
    if (is(object@resp,"lmerResp")) {
      ## result will be matrix  n x nsim :
      val <- etasim.fix + sigma * (etasim.reff +
        ## residual contribution:
        matrix(rnorm(n * nsim), ncol = nsim))
    } else if (is(object@resp,"glmResp")) {
      ## GLMM
      ## n.b. DON'T scale random-effects (???)
      	      etasim <- etasim.fix+etasim.reff
              ## FIXME:: try to avoid @call ...
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
                                  f <- factor(levels(resp)[Y+1],levels=levels(resp))
                                  split(f, rep(seq_len(nsim), each = n))
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
      tt <- terms.formula(formula(x,fixed.only=TRUE))
      attr(tt,"predvars") <- attr(attr(x@frame,"terms"),"predvars.fixed")
      tt
  }
  else attr(x@frame,"terms")
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

###----- Printing etc ----------------------------

methTitle <- function(object, dims = object@devcomp$dims)
    paste(switch(1L + dims[["GLMM"]] * 2L + dims[["NLMM"]],
                 "Linear", "Nonlinear",
                 "Generalized linear", "Generalized nonlinear"),
          "mixed model fit by",
          if(isREML(object)) "REML" else "maximum likelihood")

famlink <- function(object, resp = object@resp) {
    if(is(resp, "glmResp"))
	resp$family[c("family", "link")]
    else list(family = NULL, link = NULL)
}

.prt.family <- function(famL) {
    if (!is.null(f <- famL$family)) {
	cat(" Family:", f)
        if (!(is.null(ll <- famL$link))) cat(" (", ll, ")")
        cat("\n")
    }
}

.prt.call <- function(call) {
    if (!is.null(cc <- call$formula))
	cat("Formula:", deparse(cc),"\n")
    if (!is.null(cc <- call$data))
	cat("   Data:", deparse(cc), "\n")
    if (!is.null(cc <- call$subset))
	cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
}

getLlikAIC <- function(object, cmp = object@devcomp$cmp) {
    llik <- logLik(object)   # returns NA for a REML fit - maybe change?
    AICstats <- {
	if(isREML(object)) cmp["REML"] # *no* likelihood stats here
	else {
	    c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
	      deviance = deviance(object))
	}
    }
    list(logLik=llik, AICtab = AICstats)
}

.prt.aictab <- function(aictab, digits=4) {
    t.4 <- round(aictab, digits)
    if (length(aictab) == 1 && names(aictab) == "REML")
	cat("REML criterion at convergence:", t.4, "\n")
    else print(t.4)
}

.prt.VC <- function(varcor, digits, comp, ...) {
    cat("Random effects:\n")
    fVC <- if(missing(comp))
	formatVC(varcor, digits=digits)
    else
	formatVC(varcor, digits=digits, comp=comp)
    print(fVC, quote = FALSE, digits = digits, ...)
}

.prt.grps <- function(ngrps, nobs) {
    cat(sprintf("Number of obs: %d, groups: ", nobs))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
}

## This is modeled a bit after	print.summary.lm :
## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
			correlation = NULL, symbolic.cor = FALSE,
			signif.stars = getOption("show.signif.stars"),
			ranef.comp = c("Variance", "Std.Dev."), ...)
{
    so <- summary(x)
    cat(sprintf("%s ['%s']\n",so$methTitle, so$objClass))
    .prt.family(so)
    ## FIXME: commenting out for now, restore after release?
    ## cat("Scaled residuals:\n")
    ## print(summary(residuals(x,type="pearson",scaled=TRUE)),digits=digits)
    .prt.call(so$call); cat("\n")
    .prt.aictab(so$AICtab, 4); cat("\n")
    .prt.VC(so$varcor, digits=digits, useScale= so$useScale,
	    comp = ranef.comp, ...)
    .prt.grps(so$ngrps, nobs= so$devcomp$dims[["n"]])

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
print.merMod <- function(x, digits = max(3, getOption("digits") - 3),
                         correlation = NULL, symbolic.cor = FALSE,
                         signif.stars = getOption("show.signif.stars"),
			 ranef.comp = "Std.Dev.", ...)
{
    dims <- (devC <- x@devcomp)$dims
    methTit <- methTitle(x, dims=dims)
    cat(sprintf("%s ['%s']\n",methTit, class(x)))
    famL <- famlink(x, resp = x@resp)
    .prt.family(famL)
    .prt.call(x@call)
    useScale <- as.logical(dims[["useSc"]])

    llAIC <- getLlikAIC(x)
    .prt.aictab(llAIC$AICtab, 4)
    varcor <- VarCorr(x)
    .prt.VC(varcor, digits=digits, comp = ranef.comp, ...)
    ngrps <- sapply(x@flist, function(x) length(levels(x)))
    .prt.grps(ngrps, nobs= dims[["n"]])
    if(length(cf <- fixef(x)) >= 0) {
	cat("Fixed Effects:\n")
	print.default(format(cf, digits = digits),
		      print.gap = 2L, quote = FALSE, ...)
    } else cat("No fixed effect coefficients\n")
    invisible(x)
}
##' @exportMethod show
setMethod("show",  "merMod", function(object) print.merMod(object))

##' @S3method print summary.merMod
print.summary.merMod <- printMerenv

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

##' Construct names of individual theta/sd:cor components
##'
##' @param object a fixed model
##' @param diag.only include only diagonal elements?
##' @param old (logical) give backward-compatible results?
##' @param prefix a character vector with two elements giving the prefix
##' for diagonal (e.g. "sd") and off-diagonal (e.g. "cor") elements
##' ## @export
tnames <- function(object,diag.only=FALSE,old=TRUE,prefix=NULL) {
    if (old) {
        nc <- c(unlist(mapply(function(g,e) {
            mm <- outer(e,e,paste,sep=".")
            diag(mm) <- e
            mm <- if (diag.only) diag(mm) else mm[lower.tri(mm,diag=TRUE)]
            paste(g,mm,sep=".")
        },
        names(object@cnms),object@cnms)))
        return(nc)
    } else {
        pfun <- function(g,e) {
            mm <- outer(e,e,paste,sep=".")
            mm[] <- paste(mm,g,sep="|")
            if (!is.null(prefix)) mm[] <- paste(prefix[2],mm,sep="_")
            diag(mm) <- paste(e,g,sep="|")
            if (!is.null(prefix))  diag(mm) <- paste(prefix[1],diag(mm),sep="_")
            mm <- if (diag.only) diag(mm) else mm[lower.tri(mm,diag=TRUE)]
        }
        nc <- c(unlist(mapply(pfun,names(object@cnms),object@cnms)))
        return(nc)
    }
}

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
##' @param name a character vector specifying the name(s) of the
##' \dQuote{component}. If \code{length(name)}>1, a named list of
##' components will be returned. Possible values are:\cr
##' \describe{
##'     \item{X}{fixed-effects model matrix}
##'     \item{Z}{random-effects model matrix}
##'     \item{Zt}{transpose of random-effects model matrix.  Note that
##'              the structure of \code{Zt} has changed since \code{lme4.0};
##'              to get a backward-compatible structure, use
##'              \code{do.call(Matrix::rBind,getME(.,"Ztlist"))}}
##'     \item{Ztlist}{list of components of the transpose of the random-effects model matrix,
##'              separated by individual variance component}
##'     \item{y}{response vector}
##'     \item{mu}{conditional mean of the response}
##'     \item{u}{conditional mode of the \dQuote{spherical} random effects variable}
##'     \item{b}{conditional mode of the random effects variable}
##'     \item{Gp}{groups pointer vector.  A pointer to the beginning
##'               of each group of random effects corresponding to the
##'               random-effects terms, beginning with 0 and including
##'               a final element giving the total number of random effects}
##'     \item{Tp}{theta pointer vector.  A pointer to the beginning
##'               of the theta sub-vectors corresponding to the
##'               random-effects terms, beginning with 0 and including
##'               a final element giving the total number of random effects}
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
##'     \item{sigma}{residual standard error}
##'     \item{flist}{a list of the grouping variables (factors) involved in the random effect terms}
##'     \item{beta}{fixed-effects parameter estimates (identical to the result of \code{\link{fixef}}, but without names)}
##'     \item{theta}{random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term}
##'     \item{ST}{a list of matrices giving the relative Cholesky factors for each random effect term}
##'     \item{n_rtrms}{number of random-effects terms}
##'     \item{n_rfacs}{number of distinct random-effects grouping factors}
##'     \item{REML}{restricted maximum likelihood}
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
		  name = c("X", "Z","Zt", "Ztlist",
                  "y", "mu", "u", "b",
		  "Gp", "Tp",
		  "L", "Lambda", "Lambdat", "Lind", "A",
		  "RX", "RZX", "sigma",
                  "flist",
                  "beta", "theta", "ST",
		  "REML", "is_REML",
                  "n_rtrms", "n_rfacs", "cnms",
                  "devcomp", "offset", "lower"))
{
    if(missing(name)) stop("'name' must not be missing")
    stopifnot(is(object,"merMod"))
    if (length(name <- as.character(name))>1) {
        names(name) <- name
        return(lapply(name, getME, object=object))
    }
    name <- match.arg(name)
    rsp  <- object@resp
    PR   <- object@pp
    dc   <- object@devcomp
    cmp  <- dc $ cmp
    cnms <- object@cnms
    dims <- dc $ dims
    Tpfun <- function(cnms) {
	ltsize <- function(n) n*(n+1)/2 # lower triangle size
	cLen <- cumsum(ltsize(vapply(cnms,length, 1L)))
	setNames(c(0, cLen),
		 c("beg__", names(cnms))) ## such that diff( Tp ) is well-named
    }
    switch(name,
	   "X" = PR$X, ## ok ? - check -- use model.matrix() method instead?
	   "Z" = t(PR$Zt),
	   "Zt"= PR$Zt,
           "Ztlist" =
       {
           getInds <- function(i) {
               n <- diff(object@Gp)[i]      ## number of elements in this block
               nt <- length(cnms[[i]]) ## number of REs
               inds <- lapply(seq(nt),seq,to=n,by=nt)  ## pull out individual RE indices
               inds <- lapply(inds,function(x) x + object@Gp[i])  ## add group offset
           }
           inds <- do.call(c,lapply(seq_along(cnms),getInds))
           setNames(lapply(inds,function(i) PR$Zt[i,]),
                    tnames(object,diag.only=TRUE))
       },
           "y" = rsp$y,
           "mu"= rsp$mu,
           "u" = object@u,
           "b" = t(PR$Lambdat) %*% object@u,
	   "L"= PR$ L(),
	   "Lambda"= t(PR$ Lambdat),
	   "Lambdat"= PR$ Lambdat,
           "A" = PR$Lambdat %*% PR$Zt,
           "Lind" = PR$ Lind,
	   "RX" = structure(PR$RX(), dimnames = list(colnames(PR$X), colnames(PR$X))), ## maybe add names elsewhere?
	   "RZX" = structure(PR$RZX, dimnames = list(NULL, colnames(PR$X))), ## maybe add names elsewhere?
           "sigma" = sigma(object),
           "Gp" = object@Gp,
           "Tp" = Tpfun(cnms) ,
           "flist" = object@flist,
	   "beta" = object@beta,
           "theta"= setNames(object@theta,tnames(object)),
           "ST"= setNames(vec2STlist(object@theta,
                                    n=sapply(cnms,length)),
                         names(cnms)),
	   "REML" = dims["REML"],
	   "is_REML" = isREML(object),
           ## number of random-effects terms
	   "n_rtrms" = length(cnms),
           ## number of random-effects grouping factors
           "n_rfacs" = length(object@flist),
           "cnms" = cnms,
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
	rr@factors$correlation <-
	    if(!is.na(sigma)) as(rr, "corMatrix") else rr # (is NA anyway)
    rr
}

##' @importFrom stats vcov
##' @S3method vcov merMod
vcov.merMod <- function(object, correlation = TRUE, sigm = sigma(object), ...)
    mkVcov(sigm, unsc = object@pp$unsc(), nmsX = colnames(object@pp$X),
	   correlation=correlation, ...)

##' @importFrom stats vcov
##' @S3method vcov summary.merMod
vcov.summary.merMod <- function(object, correlation = TRUE, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of merMod object")
    object$vcov
}

##' Make variance and correlation matrices from \code{theta}
##'
##' @param sc scale factor (residual standard deviation)
##' @param cnms component names
##' @param nc numeric vector: number of terms in each RE component
##' @param theta theta vector (lower-triangle of Cholesky factors)
##' @param nms component names (FIXME: nms/cnms redundant: nms=names(cnms)?)
##' @seealso \code{\link{VarCorr}}
##' @return A matrix
##' @export
mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    if(!all(nms == names(cnms))) ## the above FIXME
	warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n",
		"Please report!", immediate.=TRUE)
    ans <- lapply(ncseq, function(i)
	      {
		  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		  Li <- diag(nrow = nc[i])
		  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		  rownames(Li) <- cnms[[i]]
		  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
		  val <- tcrossprod(sc * Li) # variance-covariance
		  stddev <- sqrt(diag(val))
		  corr <- t(val / stddev)/stddev
		  diag(corr) <- 1
		  structure(val, stddev = stddev, correlation = corr)
	      })
    if(is.character(nms)) {
	## FIXME: do we want this?  Maybe not.
	## Potential problem: the names of the elements of the VarCorr() list
	##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
	##  Subject terms, so the names are "Subject", "Subject".  The print method
	##  for VarCorrs handles this just fine, but it's a little awkward if we
	##  want to dig out elements of the VarCorr list ... ???
	if (anyDuplicated(nms))
	    nms <- make.names(nms, unique=TRUE)
	names(ans) <- nms
    }
    structure(ans, sc = sc)
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
  ## FIXME:: would like to fix nlme to add ...
  ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    if (is.null(cnms <- x@cnms))
	stop("VarCorr methods require reTrms, not just reModule")
    if(missing(sigma)) # "bug": fails via default 'sigma=sigma(x)'
	sigma <- lme4::sigma(x)  ## FIXME: do we still need lme4:: ?
    nc <- vapply(cnms, length, 1L) # no. of columns per term
    structure(mkVarCorr(sigma, cnms=cnms, nc=nc, theta = x@theta,
			nms = {fl <- x@flist; names(fl)[attr(fl, "assign")]}),
	      useSc = as.logical(x@devcomp$dims["useSc"]),
	      class = "VarCorr.merMod")
}

if(FALSE)## *NOWHERE* used _FIXME_ ??
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

##' @S3method print VarCorr.merMod
print.VarCorr.merMod <- function(x, digits = max(3, getOption("digits") - 2),
		   comp = "Std.Dev.", ...)
    print(formatVC(x, digits=digits, comp=comp), quote=FALSE, ...)

##' __NOT YET EXPORTED__
##' "format()" the 'VarCorr' matrix of the random effects -- for
##' print()ing and show()ing
##'
##' @title Format the 'VarCorr' Matrix of Random Effects
##' @param varc a \code{\link{VarCorr}} (-like) matrix with attributes.
##' @param digits the number of significant digits.
##' @param comp character vector of length one or two indicating which
##' columns out of "Variance" and "Std.Dev." should be shown in the
##' formatted output.
##' @return a character matrix of formatted VarCorr entries from \code{varc}.
formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
		     comp = "Std.Dev.")
{
    c.nms <- c("Groups", "Name", "Variance", "Std.Dev.")
    avail.c <- c.nms[-(1:2)]
    if(any(is.na(mcc <- pmatch(comp, avail.c))))
	stop("Illegal 'comp': ", comp[is.na(mcc)])
    nc <- length(colnms <- c(c.nms[1:2], (use.c <- avail.c[mcc])))
    if(length(use.c) == 0)
	stop("Must *either* show variances or standard deviations")
    useScale <- attr(varc, "useSc")
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"),
		  if(useScale) list(Residual = unname(attr(varc, "sc"))))
    reLens <- vapply(reStdDev, length, 1L)
    nr <- sum(reLens)
    reMat <- array('', c(nr, nc), list(rep.int('', nr), colnms))
    reMat[1+cumsum(reLens)-reLens, "Groups"] <- names(reLens)
    reMat[,"Name"] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
    if(any("Variance" == use.c))
    reMat[,"Variance"] <- format(unlist(reStdDev)^2, digits = digits)
    if(any("Std.Dev." == use.c))
    reMat[,"Std.Dev."] <- format(unlist(reStdDev),   digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x) {
			       x <- as(x, "matrix")
			       dig <- max(2, digits - 2) # use 'digits' !
			       cc <- format(round(x, dig), nsmall = dig)
			       cc[!lower.tri(cc)] <- ""
			       nr <- nrow(cc)
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }))[, -maxlen, drop = FALSE]
	if (nrow(corr) < nrow(reMat))
	    corr <- rbind(corr, matrix("", nrow(reMat) - nrow(corr), ncol(corr)))
	colnames(corr) <- c("Corr", rep.int("", max(0L, ncol(corr)-1L)))
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
    useSc <- as.logical(dd[["useSc"]])
    sig <- sigma(object)
    REML <- isREML(object)

    famL <- famlink(resp=resp)
    coefs <- cbind("Estimate" = fixef(object),
		   "Std. Error" = sig * sqrt(diag(object@pp$unsc())))
    if (nrow(coefs) > 0) {
	coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]), deparse.level=0)
	colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
	if (isGLMM(object))
	    coefs <- cbind(coefs, "Pr(>|z|)" =
			   2*pnorm(abs(cf3), lower.tail=FALSE))
    }

    llAIC <- getLlikAIC(object)
    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
					# use S3 class for now
    structure(list(methTitle = methTitle(object, dims=dd),
		   objClass = class(object),
		   devcomp = devC,
		   isLmer=is(resp, "lmerResp"), useScale=useSc,
		   logLik=llAIC[["logLik"]],
		   family=famL$fami, link=famL$link,
		   ngrps=sapply(object@flist, function(x) length(levels(x))),
		   coefficients=coefs, sigma=sig,
		   vcov=vcov(object, correlation=TRUE, sigm=sig),
		   varcor=varcor, # and use formatVC(.) for printing.
		   AICtab = llAIC[["AICtab"]], call=object@call
		   ), class = "summary.merMod")
}

##' @S3method summary summary.merMod
summary.summary.merMod <- function(object, varcov = TRUE, ...) {
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
	if (!is.null(pv <- attr(x, "postVar")))
        {
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
    optfun <- if (is.character(optimizer))
	tryCatch(get(optimizer), error=function(e) NULL)
    if (is.null(optfun)) stop("couldn't find optimizer function ",optimizer)
    if (!is.function(optfun)) stop("non-function specified as optimizer")
    needArgs <- c("fn","par","lower","control")
    if (any(is.na(match(needArgs, names(formals(optfun))))))
	stop("optimizer function must use (at least) formal parameters ",
	     paste(sQuote(needArgs), collapse=", "))
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
                   if (is.null(control$xt)) control$xt <- control$xst*5e-4
               })
    if (optimizer=="Nelder_Mead") control$verbose <- verbose
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
    curWarnings <- list()
    opt <- withCallingHandlers(do.call(optfun,arglist),
                               warning = function(w) {
                                   curWarnings <<- append(curWarnings,list(w$message))
                               })
    ## cat("***",unlist(tail(curWarnings,1)),"\n")
    ## FIXME: set code to warn on convergence !=0
    ## post-fit tweaking
    if (optimizer=="bobyqa") {
        opt$convergence <- opt$ierr
    }
    if (optimizer=="optimx") {
        ## optr <- lapply(opt,"[[",1)[c("par","fvalues","conv")]
        ## opt$message <- attr(opt,"details")[[1]]$message
        opt <- list(par=coef(opt)[1,],
                    fvalues=opt$value[1],
                    conv=opt$convcode[1],
                    message=attr(opt,"details")[,"message"][[1]])
    }
    if (opt$conv!=0) {
        wmsg <- paste("convergence code",opt$conv,"from",optimizer)
        if (!is.null(opt$msg)) wmsg <- paste0(wmsg,": ",opt$msg)
        warning(wmsg)
        curWarnings <<- append(curWarnings,list(wmsg))
    }
    ## store all auxiliary information
    attr(opt,"optimizer") <- optimizer
    attr(opt,"control") <- control
    attr(opt,"warnings") <- curWarnings
    opt
}
