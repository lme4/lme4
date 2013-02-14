#' @export
lformula <- function(formula, data=NULL, REML = TRUE, sparseX = FALSE,
                     control = list(), start = NULL,
                     verbose = 0L, subset, weights, na.action, offset,
                     contrasts = NULL, devFunOnly=FALSE,
                     optimizer="Nelder_Mead", ...)
{
  
  mf <- mc <- match.call()
  checkArgs("lmer",sparseX,...)
  if (!is.null(list(...)[["family"]])) {
    ## lmer(...,family=...); warning issued within checkArgs
    mc[[1]] <- as.name("glformula")
    return(eval(mc, parent.frame()) )
  }
  
  denv <- checkFormulaData(formula,data)
  mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
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
  if (any(unlist(lapply(reTrms$flist, nlevels)) >= nrow(fr)))
    stop("number of levels of each grouping factor must be ",
         "less than number of obs")
  ## fixed-effects model matrix X - remove random effects from formula:
  fixedform <- formula
  fixedform[[3]] <- if(is.null(nb <- nobars(fixedform[[3]]))) 1 else nb
  mf$formula <- fixedform
  ## re-evaluate model frame to extract predvars component
  fixedfr <- eval(mf, parent.frame())
  attr(attr(fr,"terms"),"predvars.fixed") <- attr(attr(fixedfr,"terms"),"predvars")
  X <- model.matrix(fixedform, fixedfr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
  p <- ncol(X)
  if ((rankX <- rankMatrix(X)) < p)
    stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))
  return(list(fr = fr, X = X, reTrms = reTrms, REML = REML, start = start))
}

#' @export
mklmerdevfun <- function(fr, X, reTrms, REML = TRUE, start = NULL){
  p <- ncol(X) # maybe also do rank check on X here??
  rho <- new.env(parent=parent.env(environment()))
  rho$pp <- do.call(merPredD$new, c(reTrms[c("Zt","theta","Lambdat","Lind")], n=nrow(X), list(X=X)))
  rho$resp <- mkRespMod(fr, if(REML) p else 0L) 
  # note:  REML does double duty as rank of X and a flag for using REML
  # maybe this should be mentioned in the help file for mkRespMod??
  # currently that help file says REML is logical
  devfun <- mkdevfun(rho, 0L)
  ## FIXME: should we apply start elsewhere? what about reTrms$theta?
  if (!is.null(start)) {
    rho$pp$setTheta(unlist(start))
  }
  devfun(rho$pp$theta) # one evaluation to ensure all values are set
  return(devfun) # this should pass the rho environment implicitly
}

#' @export
optimizelmer <- function(devfun, reTrms, 
                         optimizer="Nelder_Mead", verbose = 0L,
                         control = list()){
  verbose <- as.integer(verbose)
  restart <- TRUE ## FIXME; set default elsewhere?
  if (!is.null(control$restart)) {
    restart <- control$restart
    control$restart <- NULL
  }
  rho <- environment(devfun)
  opt <- optwrap(optimizer,
                 devfun, rho$pp$theta, lower=reTrms$lower, control=control,
                 adj=FALSE, verbose=verbose)
  
  if (restart) {
    ## FIXME: should we be looking at rho$pp$theta or opt$par
    ##  at this point???  in koller example (for getData(13)) we have
    ##   rho$pp$theta=0, opt$par=0.08
    if (length(bvals <- which(rho$pp$theta==reTrms$lower))>0) {
      ## *don't* use numDeriv -- cruder but fewer dependencies, no worries
      ##  about keeping to the interior of the allowed space
      theta0 <- new("numeric",rho$pp$theta) ## 'deep' copy ...
      d0 <- devfun(theta0)
      btol <- 1e-5  ## FIXME: make user-settable?
      ## FIXME: opt$fval is wrong
      bgrad <- sapply(bvals,
                      function(i) {
                        bndval <- reTrms$lower[i]
                        theta <- theta0
                        theta[i] <- bndval+btol
                        (devfun(theta)-d0)/btol
                      })
      ## what do I need to do to reset rho$pp$theta to original value???
      devfun(theta0) ## reset rho$pp$theta after tests
      ## FIXME: allow user to specify ALWAYS restart if on boundary?
      if (any(bgrad<0)) {
        if (verbose) message("some theta parameters on the boundary, restarting")
        opt <- optwrap(optimizer,
                       devfun, opt$par,
                       lower=reTrms$lower, control=control,
                       adj=FALSE, verbose=verbose)
      }
    }
    return(opt)
  }
}

#' @export
lmer2 <- function(formula, data=NULL, REML = TRUE, sparseX = FALSE,
                  control = list(), start = NULL,
                  verbose = 0L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly=FALSE,
                  optimizer="Nelder_Mead", ...){
  mc <- mcout <- match.call()
  
  # parse data and formula
  mc[[1]] <- as.name("lformula")
  lmod <- eval(mc, parent.frame(2L))
  
  # create deviance function for covariance parameters (theta)
  devfun <- do.call(mklmerdevfun, lmod)
  
  # optimze deviance function over covariance parameters
  opt <- optimizelmer(devfun, lmod$reTrms, ...)
  
  # prepare output
  mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr, mcout)
}

## TODO:  remove any arguments that aren't actually used by glformula (same for lformula)
#' @export
glformula <- function(formula, data=NULL, family = gaussian, sparseX = FALSE,
                      control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                      compDev = TRUE, subset, weights, na.action, offset,
                      contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
                      tolPwrss = 1e-7, optimizer=c("bobyqa","Nelder_Mead"), ...){
  ## FIXME: probably don't need devFunOnly
  ## FIXME: does start= do anything? test & fix
  verbose <- as.integer(verbose)
  mf <- mc <- match.call()
  # extract family, call lmer for gaussian
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    mc[[1]] <- as.name("lformula")
    mc["family"] <- NULL            # to avoid an infinite loop
    return(eval(mc, parent.frame()))
  }
  if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
    stop('"quasi" families cannot be used in glmer')
  
  checkArgs("glmer",sparseX,...)
  
  stopifnot(length(nAGQ <- as.integer(nAGQ)) == 1L,
            nAGQ >= 0L,
            nAGQ <= 25L)
  
  denv <- checkFormulaData(formula,data)
  mc$formula <- formula <- as.formula(formula,env=denv)    ## substitute evaluated version
  
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
  
  if ((rankX <- rankMatrix(X)) < p)
    stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))
  
  return(list(fr = fr, X = X, reTrms = reTrms, family = family))
}

#' @export
mkglmerdevfun <- function(fr, X, reTrms, family,
                          nAGQ = 1L, verbose = 0L, tolPwrss = 1e-7, compDev = TRUE){
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
  #rho$GQmat       <- GHrule(nAGQ) # experimental!! SCW
  devfun <- mkdevfun(rho, 0L) # changing 0L to 1L is also experimental!! SCW
  #if (devFunOnly && !nAGQ) return(devfun)
  return(devfun) # this should pass the rho environment implicitly
}

#' @export
optimizeglmer <- function(devfun, optimizer="bobyqa", 
                          nAGQ = 1L, verbose = 0L, control = list()){
  rho <- environment(devfun)
  opt <- optwrap(optimizer, devfun,rho$pp$theta, rho$lower,
                 control=control, adj=FALSE, verbose=verbose)
  
  # not sure where next two lines should be placed??  maybe mkglmerdevfun2??
  # put here b/c i want to avoid an if structure in mkglmerdevfun2
  rho$control <- attr(opt,"control")
  rho$nAGQ <- nAGQ
  
  return(opt)
}

## only do this function if nAGQ > 0L
#' @export
updateglmerdevfun <- function(devfun, reTrms, nAGQ = 1L){
  rho <- environment(devfun)
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
  devfun <- mkdevfun(rho, nAGQ)  # does this attach rho to devfun??
  return(devfun)
}

#' @export
reoptimizeglmer <- function(devfun, optimizer="Nelder_Mead", 
                            verbose = 0L, control = list()){
  rho <- environment(devfun)
  opt <- optwrap(optimizer,devfun,c(rho$pp$theta, rho$pp$delb),
                 rho$lower, control=control,
                 adj=TRUE, verbose=verbose)
  rho$resp$setOffset(rho$baseOffset)
  return(opt)
}

#' @export
glmer2 <- function(formula, data=NULL, family = gaussian, sparseX = FALSE,
                   control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                   compDev = TRUE, subset, weights, na.action, offset,
                   contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
                   tolPwrss = 1e-7, optimizer=c("bobyqa","Nelder_Mead"), ...){
  mc <- mcout <- match.call()
  
  # parse the formula and data
  mc[[1]] <- as.name("glformula")
  glmod <- eval(mc, parent.frame(2L))
  
  # create deviance function for covariance parameters (theta)
  devfun <- do.call(mkglmerdevfun, c(glmod, list(compDev = compDev)))
  
  # optimize deviance function over covariance parameters
  opt <- optimizeglmer(devfun, optimizer = optimizer[[1]], ...)
  
  if(nAGQ > 0L){
    # update deviance function to include fixed effects as inputs
    devfun <- updateglmerdevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
    # reoptimize deviance function over covariance parameters and fixed effects
    opt <- reoptimizeglmer(devfun, optimizer = optimizer[[2]], 
                           verbose = verbose, control = control)
  }
  
  # prepare output
  mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr, mcout)
}
