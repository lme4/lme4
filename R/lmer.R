# lmer, glmer and nlmer plus methods and utilities

if (FALSE) {
### FIXME: Move this function to the stats package
rWishart <- function(n, df, invScal)
### Random sample from a Wishart distribution
    .Call(lme4_rWishart, n, df, invScal)
}

### Utilities for parsing the mixed model formula

findbars <- function(term)
### Return the pairs of expressions that separated by vertical bars
{
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

nobars <- function(term)
### Return the formula omitting the pairs of expressions that are
### separated by vertical bars
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

subbars <- function(term)
### Substitute the '+' function for the '|' function
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

subnms <- function(term, nlist)
### Substitute any names from nlist in term with 1
{
    if (!is.language(term)) return(term)
    if (is.name(term)) {
        if (any(unlist(lapply(nlist, get("=="), term)))) return(1)
        return(term)
    }
    stopifnot(length(term) >= 2)
    for (j in 2:length(term)) term[[j]] <- subnms(term[[j]], nlist)
    term
}

slashTerms <- function(x)
### Return the list of '/'-separated terms in an expression that
### contains slashes
{
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

makeInteraction <- function(x)
### from a list of length 2 return recursive interaction terms
{
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
}

expandSlash <- function(bb)
### expand any slashes in the grouping factors returned by findbars
{
    if (!is.list(bb)) return(expandSlash(list(bb)))
    ## I really do mean lapply(unlist(... - unlist returns a
    ## flattened list in this case
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
            return(lapply(unlist(makeInteraction(trms)),
                          function(trm) substitute(foo|bar,
                                                   list(foo = x[[2]],
                                                        bar = trm))))
        x
    }))
}

### Utilities used in lmer, glmer and nlmer

createCm <- function(A, s)
### Create the nonzero pattern for the sparse matrix Cm from A.
### ncol(A) is s * ncol(Cm).  The s groups of ncol(Cm) consecutive
### columns in A are overlaid to produce Cm.
{
    stopifnot(is(A, "dgCMatrix"))
    s <- as.integer(s)[1]
    if (s == 1L) return(A)
    if ((nc <- ncol(A)) %% s)
        stop(gettextf("ncol(A) = %d is not a multiple of s = %d",
                      nc, s))
    ncC <- as.integer(nc / s)
    TA <- as(A, "TsparseMatrix")
    as(new("dgTMatrix", Dim = c(nrow(A), ncC),
           i = TA@i, j = as.integer(TA@j %% ncC), x = TA@x),
       "CsparseMatrix")
}

### FIXME: somehow the environment of the mf formula does not have
### .globalEnv in its parent list.  example(Mmmec, package = "mlmRev")
### used to have a formula of ~ offset(log(expected)) + ... and the
### offset function was not found in eval(mf, parent.frame(2))
lmerFrames <- function(mc, formula, contrasts, vnms = character(0))
### Create the model frame, X, Y, wts, offset and terms

### mc - matched call of calling function
### formula - two-sided formula
### contrasts - contrasts argument
### vnms - names of variables to be included in the model frame
{
    mf <- mc
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]

    ## The model formula for evaluation of the model frame.  It looks
    ## like a linear model formula but includes any random effects
    ## terms and any names of parameters used in a nonlinear mixed model.
    frame.form <- subbars(formula)      # substitute `+' for `|'
    if (length(vnms) > 0)               # add the variables names for nlmer
        frame.form[[3]] <-
            substitute(foo + bar,
                       list(foo = parse(text = paste(vnms, collapse = ' + '))[[1]],
                            bar = frame.form[[3]]))

    ## The model formula for the fixed-effects terms only.
    fixed.form <- nobars(formula)       # remove any terms with `|'
    if (!inherits(fixed.form, "formula"))
      ## RHS is empty - use `y ~ 1'
      fixed.form <- as.formula(substitute(foo ~ 1, list(foo = fixed.form)))

    ## attach the correct environment
    environment(fixed.form) <- environment(frame.form) <- environment(formula)

    ## evaluate a model frame
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fe <- mf                            # save a copy of the call
    mf <- eval(mf, parent.frame(2))

    ## evaluate the terms for the fixed-effects only (used in anova)
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2)) # allow model.frame to update them

    ## response vector
    Y <- model.response(mf, "any")
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
    mt <- attr(fe, "terms")

    ## Extract X checking for a null model. This check shouldn't be
    ## needed because an empty formula is changed to ~ 1 but it can't hurt.
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts) else matrix(,NROW(Y),0)
    storage.mode(X) <- "double"      # when ncol(X) == 0, X is logical
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL

    ## Extract the weights and offset.  For S4 classes we want the
    ## `not used' condition to be numeric(0) instead of NULL
    wts <- model.weights(mf); if (is.null(wts)) wts <- numeric(0)
    off <- model.offset(mf); if (is.null(off)) off <- numeric(0)

    ## check weights and offset
    if (any(wts <= 0))
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if(length(off) && length(off) != NROW(Y))
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), NROW(Y)))

    ## remove the terms attribute from mf
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), mf = mf, fixef = fixef)
}

##' Is f1 nested within f2?
##'
##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @param f1 factor 1
##' @param f2 factor 2

##' @return TRUE if factor 1 is nested within factor 2

isNested <- function(f1, f2)
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix",
                 i = as.integer(f2) - 1L,
                 j = as.integer(f1) - 1L,
                 Dim = c(length(levels(f2)),
                 length(levels(f1)))),
             "CsparseMatrix")
    all(diff(sm@p) < 2)
}

##' dimsNames and devNames are in the package's namespace rather than
##' in the function lmerFactorList because the function sparseRasch
##' needs to access them.

dimsNames <- c("nt", "n", "p", "q", "s", "np", "LMM", "REML",
               "fTyp", "lTyp", "vTyp", "nest", "useSc", "nAGQ",
               "verb", "mxit", "mxfn", "cvg")
dimsDefault <- list(s = 1L,             # identity mechanistic model
                    mxit= 300L,         # maximum number of iterations
                    mxfn= 900L, # maximum number of function evaluations
                    verb= 0L,           # no verbose output
                    np= 0L,             # number of parameters in ST
                    LMM= 0L,            # not a linear mixed model
                    REML= 0L,         # glmer and nlmer don't use REML
                    fTyp= 2L,           # default family is "gaussian"
                    lTyp= 5L,           # default link is "identity"
                    vTyp= 1L, # default variance function is "constant"
                    useSc= 1L, # default is to use the scale parameter
                    nAGQ= 1L,                  # default is Laplace
                    cvg = 0L)                  # no optimization yet attempted

devNames <- c("ML", "REML", "ldL2", "ldRX2", "sigmaML",
              "sigmaREML", "pwrss", "disc", "usqr", "wrss",
              "dev", "llik", "NULLdev")


##' Create model matrices from r.e. terms.
##'
##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##'
##' @param formula model formula
##' @param fr: list with '$mf': model frame; '$X': .. matrix
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices
##'
##' @return a list with components named \code{"trms"}, \code{"fl"}
##'        and \code{"dims"}
lmerFactorList <- function(formula, fr, rmInt, drop)
{
    mf <- fr$mf
    ## record dimensions and algorithm settings

    ## create factor list for the random effects
    bars <- expandSlash(findbars(formula[[3]]))
    if (!length(bars)) stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars,
                 function(x)
             {
                 ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                                       list(fac = x[[3]])), mf)
                 im <- as(ff, "sparseMatrix") # transpose of indicators
		 ## Could well be that we should rather check earlier .. :
		 if(!isTRUE(validObject(im, test=TRUE)))
		     stop("invalid conditioning factor in random effect: ", format(x[[3]]))

                 mm <- model.matrix(eval(substitute(~ expr, # model matrix
                                                    list(expr = x[[2]]))),
                                    mf)
                 if (rmInt) {
                     if (is.na(icol <- match("(Intercept)", colnames(mm)))) break
                     if (ncol(mm) < 2)
                         stop("lhs of a random-effects term cannot be an intercept only")
                     mm <- mm[ , -icol , drop = FALSE]
                 }
                 ans <- list(f = ff,
                             A = do.call(rBind,
                             lapply(seq_len(ncol(mm)), function(j) im)),
                             Zt = do.call(rBind,
                             lapply(seq_len(ncol(mm)),
                                    function(j) {im@x <- mm[,j]; im})),
                             ST = matrix(0, ncol(mm), ncol(mm),
                             dimnames = list(colnames(mm), colnames(mm))))
                 if (drop) {
                     ## This is only used for nlmer models.
                     ## Need to do something more complicated for A
                     ## here.  Essentially you need to create a copy
                     ## of im for each column of mm, im@x <- mm[,j],
                     ## create the appropriate number of copies,
                     ## prepend matrices of zeros, then rBind and drop0.
                     ans$A@x <- rep(0, length(ans$A@x))
                     ans$Zt <- drop0(ans$Zt)
                 }
                 ans
             })
    dd <-
        VecFromNames(dimsNames, "integer",
                     c(list(n = nrow(mf), p = ncol(fr$X), nt = length(fl),
                            q = sum(sapply(fl, function(el) nrow(el$Zt)))),
                       dimsDefault))
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true
    nlev <- sapply(fl, function(el) length(levels(el$f)))
    ## determine the number of random effects at this point
    if (any(diff(nlev)) > 0) fl <- fl[rev(order(nlev))]
    ## separate the terms from the factor list
    trms <- lapply(fl, "[", -1)
    names(trms) <- NULL
    fl <- lapply(fl, "[[", "f")
    attr(fl, "assign") <- seq_along(fl)
    ## check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        fl <- fl[match(ufn, fnms)]
        attr(fl, "assign") <- match(fnms, ufn)
    }
    names(fl) <- ufn
    ## check for nesting of factors
    dd["nest"] <- all(sapply(seq_along(fl)[-1],
                             function(i) isNested(fl[[i-1]], fl[[i]])))

    list(trms = trms, fl = fl, dims = dd)
}

checkSTform <- function(ST, STnew)
### Check that the 'STnew' argument matches the form of ST.
{
    stopifnot(is.list(STnew), length(STnew) == length(ST),
              all.equal(names(ST), names(STnew)))
    lapply(seq_along(STnew), function (i)
           stopifnot(class(STnew[[i]]) == class(ST[[i]]),
                     all.equal(dim(STnew[[i]]), dim(ST[[i]]))))
    all(unlist(lapply(STnew, function(m) all(diag(m) > 0))))
}

lmerControl <- function(msVerbose = getOption("verbose"),
                        maxIter = 300L, maxFN = 900L)
### Control parameters for lmer, glmer and nlmer
{
    stopifnot(maxIter >= 0, maxFN >= 0)
    list(
         maxIter = as.integer(maxIter),
         maxFN = as.integer(maxFN),
	 msVerbose = as.integer(msVerbose))# "integer" on purpose
}

VecFromNames <- function(nms, mode = "numeric", defaults = list())
### Generate a named vector of the given mode
{
    ans <- vector(mode = mode, length = length(nms))
    names(ans) <- nms
    ans[] <- NA
    if ((nd <- length(defaults <- as.list(defaults))) > 0) {
        if (length(dnms <- names(defaults)) < nd)
            stop("defaults must be a named list")
        stopifnot(all(dnms %in% nms))
        ans[dnms] <- as(unlist(defaults), mode)
    }
    ans
}

mkZt <- function(FL, start, s = 1L)
### Create the standard versions of flist, Zt, Gp, ST, A, Cm,
### Cx, and L. Update dd.
{
    dd <- FL$dims
    fl <- FL$fl
    asgn <- attr(fl, "assign")
    trms <- FL$trms
    ST <- lapply(trms, `[[`, "ST")
    Ztl <- lapply(trms, `[[`, "Zt")
    Zt <- do.call(rBind, Ztl)
    Zt@Dimnames <- vector("list", 2)
    Gp <- unname(c(0L, cumsum(sapply(Ztl, nrow))))
    .Call(mer_ST_initialize, ST, Gp, Zt)
    A <- do.call(rBind, lapply(trms, `[[`, "A"))
    rm(Ztl, FL)                         # because they could be large
    nc <- sapply(ST, ncol)         # of columns in els of ST
    Cm <- createCm(A, s)
    L <- .Call(mer_create_L, Cm)
    if (s < 2) Cm <- new("dgCMatrix")
    if (!is.null(start) && checkSTform(ST, start)) ST <- start

    nvc <- sapply(nc, function (qi) (qi * (qi + 1))/2) # no. of var. comp.
### FIXME: Check number of variance components versus number of
### levels in the factor for each term. Warn or stop as appropriate

    dd["np"] <- as.integer(sum(nvc))    # number of parameters in optimization
    dev <- VecFromNames(devNames, "numeric")
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn

    list(Gp = Gp, ST = ST, A = A, Cm = Cm, L = L, Zt = Zt,
         dd = dd, dev = dev, flist = fl)
}

famNms <- c("binomial", "gaussian", "Gamma", "inverse.gaussian",
            "poisson")
linkNms <- c("logit", "probit", "cauchit", "cloglog", "identity",
	     "log", "sqrt", "1/mu^2", "inverse")
varNms <- c("constant", "mu(1-mu)", "mu", "mu^2", "mu^3")

famType <- function(family)
{
    if (!(fTyp <- match(family$family, famNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    if (!(lTyp <- match(family$link, linkNms, nomatch = 0)))
        stop(gettextf("unknown link: %s",
                      sQuote(family$link), domain = "R-lme4"))
    vNam <- switch(fTyp,
                   "mu(1-mu)",          # binomial
                   "constant",          # gaussian
                   "mu^2",              # Gamma
                   "mu^3",              # inverse.gaussian
                   "mu")                # poisson
    if (!(vTyp <- match(vNam, varNms, nomatch = 0)))
        stop(gettextf("unknown GLM family: %s",
                      sQuote(family$family), domain = "R-lme4"))
    c(fTyp = fTyp, lTyp = lTyp, vTyp = vTyp)
}

convergenceMessage <- function(cvg)
### Create the convergence message
{
    msg <- switch(as.character(cvg),
                  "3" = "X-convergence (3)",
                  "4" = "relative convergence (4)",
                  "5" = "both X-convergence and relative convergence (5)",
                  "6" = "absolute function convergence (6)",

                  "7" = "singular convergence (7)",
                  "8" = "false convergence (8)",
                  "9" = "function evaluation limit reached without convergence (9)",
                  "10" = "iteration limit reached without convergence (9)",
                  "14" = "storage has been allocated (?) (14)",

                  "15" = "LIV too small (15)",
                  "16" = "LV too small (16)",
                  "63" = "fn cannot be computed at initial par (63)",
                  "65" = "gr cannot be computed at initial par (65)")
    if (is.null(msg))
        msg <- paste("See PORT documentation.  Code (", cvg, ")", sep = "")
    msg
}

mer_finalize <- function(ans)
{
    .Call(mer_optimize, ans)
    if (ans@dims[["cvg"]] > 6) warning(convergenceMessage(ans@dims[["cvg"]]))
    .Call(mer_update_ranef, ans)
    .Call(mer_update_mu, ans)
    ans
}

## Modifications to lmer often involve modifying model matrices before
## creating and optimizing the mer object.  Everything past the model
## matrices is encapsulated in this function
lmer_finalize <- function(fr, FL, start, REML, verbose)
{
    Y <- as.double(fr$Y)
    if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
        start <- list(ST = start)
    if (is.numeric(start)) start <- list(STpars = start)
    dm <- mkZt(FL, start[["ST"]])
### This checks that the number of levels in a grouping factor < n
### Only need to check the first factor because it is the one with
### the most levels.
    if (!(length(levels(dm$flist[[1]])) < length(Y)))
        stop(paste("Number of levels of a grouping factor for the random effects",
                   "must be less than the number of observations", sep = "\n"))

    dm$dd["REML"] <- as.logical(REML)
    dm$dd["verb"] <- as.integer(verbose)
    swts <- sqrt(unname(fr$wts))
    p <- dm$dd[["p"]]
    n <- length(Y)

    ans <- new(Class = "mer",
               env = new.env(),
               nlmodel = (~I(x))[[2]],
               frame = fr$mf,
               call = call("foo"),      # later overwritten
               flist = dm$flist,
               X = fr$X,
               Zt = dm$Zt,
               pWt = unname(fr$wts),
               offset = unname(fr$off),
### FIXME: Should y retain its names? As it stands any row names in the
### frame are dropped.  Really?  Are they part of the frame slot (if not
### reduced to 0 rows)?
               y = unname(Y),
               Gp = unname(dm$Gp),
               dims = dm$dd,
               ST = dm$ST,
               A = dm$A,
               Cm = dm$Cm,
	       Cx = if (length(swts)) (dm$A)@x else numeric(0),
               L = dm$L,
               deviance = dm$dev,
               fixef = fr$fixef,
               ranef = numeric(dm$dd[["q"]]),
               u = numeric(dm$dd[["q"]]),
               eta = numeric(n),
               mu = numeric(n),
               resid = numeric(n),
               sqrtrWt = swts,
               sqrtXWt = as.matrix(swts),
               RZX = matrix(0, dm$dd[["q"]], p),
               RX = matrix(0, p, p))
    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
        STp <- .Call(mer_ST_getPars, ans)
        if (length(STp) == length(stp))
            .Call(mer_ST_setPars, ans, stp)
    }
    mer_finalize(ans)
}

glmer_finalize <- function(fr, FL, glmFit, start, nAGQ, verbose)
{
    if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
        start <- list(ST = start)
    if (is.numeric(start)) start <- list(STpars = start)
    dm <- mkZt(FL, start[["ST"]])
    ft <- famType(glmFit$family)
    dm$dd[names(ft)] <- ft
    useSc <- as.integer(!(famNms[dm$dd[["fTyp"]] ] %in%
			  c("binomial", "poisson")))
    dm$dd[["useSc"]] <- useSc
    ## Only need to check the first factor because it is the one with
    ## the most levels.
    M1 <- length(levels(dm$flist[[1]]))
    n <- ncol(dm$Zt)
    if (M1 >= n) {
	msg1 <- "Number of levels of a grouping factor for the random effects\n"
	msg3 <- "n, the number of observations"
	if (useSc)
	    stop(msg1, "must be less than ", msg3)
	else if (M1 == n)
	    message(msg1, "is *equal* to ", msg3)
    }
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
    if (nAGQ %% 2 == 0) nAGQ <- nAGQ + 1L # reset nAGQ to be an odd number
    dm$dd["nAGQ"] <- as.integer(nAGQ)
    AGQlist <- .Call(lme4_ghq, nAGQ)
    y <- unname(as.double(glmFit$y))
    ##    dimnames(fr$X) <- NULL
    p <- dm$dd[["p"]]
    dm$dd["verb"] <- as.integer(verbose)
    fixef <- fr$fixef
    fixef[] <- coef(glmFit)
    if (!is.null(ff <- start$fixef) && is.numeric(ff) &&
        length(ff) == length(fixef)) fixef <- ff

    ans <- new(Class = "mer",
               env = new.env(),
               nlmodel = (~I(x))[[2]],
               frame = fr$mf,
               call = call("foo"),      # later overwritten
               flist = dm$flist,
               Zt = dm$Zt, X = fr$X, y = y,
               pWt = unname(glmFit$prior.weights),
               offset = unname(fr$off),
               Gp = unname(dm$Gp),
               dims = dm$dd, ST = dm$ST, A = dm$A,
               Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L,
               deviance = dm$dev,
               fixef = fixef,
	       ranef = numeric(dm$dd[["q"]]),
	       u = numeric(dm$dd[["q"]]),
               eta = unname(glmFit$linear.predictors),
               mu = unname(glmFit$fitted.values),
	       muEta = numeric(dm$dd[["n"]]),
	       var = numeric(dm$dd[["n"]]),
               resid = unname(glmFit$residuals),
	       sqrtXWt = as.matrix(numeric(dm$dd[["n"]])),
	       sqrtrWt = numeric(dm$dd[["n"]]),
	       RZX = matrix(0, dm$dd[["q"]], p),
               RX = matrix(0, p, p),
	       ghx = AGQlist[[1]],
	       ghw = AGQlist[[2]])
    if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
        STp <- .Call(mer_ST_getPars, ans)
        if (length(STp) == length(stp))
            .Call(mer_ST_setPars, ans, stp)
    }
    mer_finalize(ans)
    ans
}

### The main event
lmer <-
    function(formula, data, family = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
### Linear Mixed-Effects in R
{
    mc <- match.call()
    if (!is.null(family)) {             # call glmer
        mc[[1]] <- as.name("glmer")
        return(eval.parent(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lmerFactorList(formula, fr, 0L, 0L) # flist, Zt, dims
    largs <- list(...)
    if (!is.null(method <- largs$method)) {
        warning(paste("Argument", sQuote("method"),
                      "is deprecated.  Use", sQuote("REML"),
                      "instead"))
        REML <- match.arg(method, c("REML", "ML")) == "REML"
        largs <- largs[names(largs) != "method"]
    }
    if(length(largs))
	warning("the following '...' arguments have  *not* been used: ",
		sub("^list", "", deparse(largs, control=NULL)))
### FIXME: issue a warning if the control argument has an msVerbose component
    cv <- do.call(lmerControl, control)
    if (missing(verbose)) verbose <- cv$msVerbose
    FL$dims["mxit"] <- cv$maxIter
    FL$dims["mxfn"] <- cv$maxFN
    ans <- list(fr = fr, FL = FL, start = start, REML = REML, verbose = verbose)
    if (doFit) {
        ans <- do.call(lmer_finalize, ans)
        ans@call <- mc
    }
    ans
}

## for backward compatibility
lmer2 <-
    function(formula, data, family = NULL, REML = TRUE,
             control = list(), start = NULL, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    .Deprecated("lmer")
    mc <- match.call()
    mc[[1]] <- as.name("lmer")
    eval.parent(mc)
}

glmer <-
function(formula, data, family = gaussian, start = NULL,
         verbose = FALSE, nAGQ = 1, doFit = TRUE, subset, weights,
         na.action, offset, contrasts = NULL, model = TRUE,
         control = list(), ...)
### Fit a generalized linear mixed model
{
    mc <- match.call()
    ## Evaluate and check the family [[hmm.. have  famType() for that ...]]
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if(!is.list(family) || is.null(family$family))
	stop(gettextf("family '%s' not recognized", deparse(substitute(family)),
		      domain = "R-lme4"))
    if(family$family == "gaussian" && family$link == "identity") {
        mc[[1]] <- as.name("lmer")      # use lmer not glmer
        mc$family <- NULL
        return(eval.parent(mc))
    }
    if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
        stop('"quasi" families cannot be used in glmer')
    stopifnot(length(formula <- as.formula(formula)) == 3)

    ## Check for method argument which is no longer used
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
    }

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    offset <- wts <- NULL
    if (length(fr$wts)) wts <- fr$wts
    if (length(fr$off)) offset <- fr$off
    glmFit <- glm.fit(fr$X, fr$Y, weights = wts, # glm on fixed effects
                      offset = offset, family = family,
                      intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
    FL <- lmerFactorList(formula, fr, 0L, 0L) # flist, Zt
### FIXME: issue a warning if the control argument has an msVerbose component
    cv <- do.call(lmerControl, control)
    if (missing(verbose)) verbose <- cv$msVerbose
### FIXME: issue a warning if the model argument is FALSE.  It is ignored.
    FL$dims["mxit"] <- cv$maxIter
    FL$dims["mxfn"] <- cv$maxFN
    ans <- list(fr = fr, FL = FL, glmFit = glmFit, start = start,
                nAGQ = nAGQ, verbose = verbose)
    if (doFit) {
        ans <- do.call(glmer_finalize, ans)
        ans@call <- mc
    }
    ans
}

nlmer <- function(formula, data, start = NULL, verbose = FALSE,
                  nAGQ = 1, doFit = TRUE, subset, weights, na.action,
                  contrasts = NULL, model = TRUE, control = list(), ...)
### Fit a nonlinear mixed-effects model
{
    mc <- match.call()
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
    if (length(nlform) < 3)
        stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
    if (is.numeric(start)) start <- list(fixef = start)
    s <- length(pnames <- names(start$fixef))
    stopifnot(length(start$fixef) > 0, s > 0,
              inherits(data, "data.frame"), nrow(data) > 1)
### FIXME: Allow for the data argument to be missing.  What should the
### default be?
    if (any(pnames %in% names(data)))
        stop("parameter names must be distinct from names of the variables in data")
    anms <- all.vars(nlmod)
    if (!all(pnames %in% anms))
        stop("not all parameter names are used in the nonlinear model expression")

    if (!length(vnms <- setdiff(anms, pnames)))
        stop("there are no variables used in the nonlinear model expression")
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L

    ## create a frame in which to evaluate the factor list
    fr <- lmerFrames(mc,
                     eval(substitute(foo ~ bar,
                                     list(foo = nlform[[2]],
                                          bar = subnms(formula[[3]],
                                          lapply(pnames, as.name))))),
                     contrasts, vnms)
    mf <- fr$mf
    env <- new.env()
    lapply(names(mf), function(nm) assign(nm, env = env, mf[[nm]]))
    n <- nrow(mf)
    lapply(pnames,
           function(nm) assign(nm, env = env, rep(start$fixef[[nm]],
                                   length.out = n)))

    n <- nrow(mf)
    mf <- mf[rep(seq_len(n), s), ]
    row.names(mf) <- NULL
    ss <- rep.int(n, s)
    for (nm in pnames)
        mf[[nm]] <- rep.int(as.numeric(nm == pnames), ss)
    fr$mf <- mf
                                        # factor list and model matrices
    FL <- lmerFactorList(substitute(foo ~ bar, list(foo = nlform[[2]],
                                                    bar = formula[[3]])),
                         fr, TRUE, TRUE)
    X <- as.matrix(mf[,pnames])
    rownames(X) <- NULL
    xnms <- colnames(fr$X)
    if (!is.na(icol <- match("(Intercept)",xnms))) xnms <- xnms[-icol]
### FIXME: The only times there would be additional columns in the
### fixed effects would be as interactions with parameter names and
### they must be constructed differently
#    if (length(xnms) > 0)
#        Xt <- cbind(Xt, fr$X[rep.int(seq_len(n), s), xnms, drop = FALSE])
    dm <- mkZt(FL, start$STpars, s)
    cv <- do.call("lmerControl", control)
    if (missing(verbose)) verbose <- cv$msVerbose
    dm$dd["verb"] <- as.integer(verbose)
    p <- dm$dd["p"] <- length(start$fixef)
### FIXME: It is better to have lmerFactorList take the value of s
    dm$dd["n"] <- n
    dm$dd["s"] <- s
    if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
    if (nAGQ %% 2 == 0) nAGQ <- nAGQ + 1L      # reset nAGQ to be an odd number
    dm$dd["nAGQ"] <- nAGQ
    AGQlist = .Call(lme4_ghq, nAGQ)

    ans <- new(Class = "mer",
               env = env,
               nlmodel = nlmod,
               frame = fr$mf,
               call = mc,
               flist = dm$flist,
               X = X,
               Zt = dm$Zt,
               pWt = unname(sqrt(fr$wts)),
               offset = unname(fr$off),
               y = unname(as.double(fr$Y)),
               Gp = unname(dm$Gp),
               dims = dm$dd,
               ## slots that change during the iterations
               ST = dm$ST,
               V = matrix(0, n, s, dimnames = list(NULL, pnames)),
               A = dm$A,
               Cm = dm$Cm,
               L = dm$L,
               deviance = dm$dev,
               fixef = start$fixef,
	       ranef = numeric(dm$dd[["q"]]),
	       u = numeric(dm$dd[["q"]]),
               eta = numeric(n),
               mu = numeric(n),
               resid = numeric(n),
               sqrtXWt = matrix(0, n, s, dimnames = list(NULL, pnames)),
               sqrtrWt = unname(sqrt(fr$wts)),
               RZX = matrix(0, dm$dd[["q"]], p),
               RX = matrix(0, p, p),
	       ghx = AGQlist[[1]],
	       ghw = AGQlist[[2]]
               )
    .Call(mer_update_mu, ans)
### Add a check that the parameter names match the column names of gradient
    mer_finalize(ans)
}

#### Extractors specific to mixed-effects models

coef.mer <- function(object, ...)
{
    if (length(list(...)))
        warning(paste('arguments named "',
                      paste(names(list(...)), collapse = ", "),
                      '" ignored', sep = ''))
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    val <- lapply(ref, function(x) fef[rep(1, nrow(x)),,drop = FALSE])
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
}

setMethod("coef", signature(object = "mer"), coef.mer)
setMethod("coef", signature(object = "summary.mer"),
          function(object, ...) object@coefs)
## questionable whether this should be added
#setMethod("coefficients", signature(object = "mer"), coef.mer)


setAs("mer", "dtCMatrix", function(from)
### Extract the L matrix
      as(from@L, "sparseMatrix"))

setMethod("fixef", signature(object = "mer"),
          function(object, ...)
### Extract the fixed effects
          object@fixef)

##' Create a list of lists from multiple parallel lists

##' @param A a list
##' @param ... other, parallel lists

##' @return a list of lists

plist <- function(A, ...)
{
    dots <- list(...)
    stopifnot(is.list(A), all(sapply(dots, is.list)),
              all(sapply(dots, length) == length(A)))
    dots <- c(list(A), dots)
    ans <- A
    for (i in seq_along(A)) ans[[i]] <- lapply(dots, "[[", i)
    ans
}

##' Extract the random effects.
##'
##' Extract the conditional modes, which for a linear mixed model are
##' also the conditional means, of the random effects, given the
##' observed responses.  These also depend on the model parameters.
##'
##' @param object an object that inherits from the \code{\linkS4class{mer}} class
##' @param postVar logical scalar - should the posterior variance be returned
##' @param drop logical scalar - drop dimensions of single extent
##' @param whichel - vector of names of factors for which to return results

##' @return a named list of arrays or vectors, aligned to the factor list

setMethod("ranef", signature(object = "mer"),
          function(object, postVar = FALSE, drop = FALSE, whichel = names(wt), ...)
      {
          rr <- object@ranef
          ## nt is the number of terms, cn is the list of column names
          nt <- length(cn <- lapply(object@ST, colnames))
          lterm <- lapply(plist(reinds(object@Gp), cn),
                          function(el) {
                              cni <- el[[2]]
                              matrix(rr[ el[[1]] ], nc = length(cni),
                                     dimnames = list(NULL, cni))
                          })
          wt <- whichterms(object)
          ans <- lapply(plist(wt, object@flist),
                        function(el) {
                            ans <- do.call(cbind, lterm[ el[[1]] ])
                            rownames(ans) <- levels(el[[2]])
                            data.frame(ans, check.names = FALSE)
                        })
          ## Process whichel
          stopifnot(is(whichel, "character"))
          whchL <- names(wt) %in% whichel
          ans <- ans[whchL]

          if (postVar) {
              pV <- .Call(mer_postVar, object, whchL)
              for (i in seq_along(ans))
                  attr(ans[[i]], "postVar") <- pV[[i]]
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
          ans
      })

print.ranef.mer <- function(x, ...) print(unclass(x), ...)
print.coef.mer <- function(x, ...) print(unclass(x), ...)

setMethod("sigma", signature(object = "mer"),
          function (object, ...) {
              dd <- object@dims
	      if(!dd[["useSc"]]) return(1)
	      object@deviance[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]]
          })

setMethod("VarCorr", signature(x = "mer"),
	  function(x, ...)
### Create the VarCorr object of variances and covariances
      {
          sc <- sigma(x)
	  ans <- lapply(cc <- .Call(mer_ST_chol, x),
                        function(ch) {
                            val <- crossprod(sc * ch) # variance-covariance
                            stddev <- sqrt(diag(val))
                            correl <- t(val / stddev)/stddev
                            diag(correl) <- 1
                            attr(val, "stddev") <- stddev
                            attr(val, "correlation") <- correl
                            val
                        })
          fl <- x@flist
          names(ans) <- names(fl)[attr(fl, "assign")]
          attr(ans, "sc") <- if (x@dims[["useSc"]]) sc else NA
          ans
      })

#### Methods for standard extractors for fitted models

setMethod("anova", signature(object = "mer"),
	  function(object, ...)
      {
	  mCall <- match.call(expand.dots = TRUE)
	  dots <- list(...)
	  modp <- if (length(dots))
	      sapply(dots, is, "mer") | sapply(dots, is, "lm") else logical(0)
	  if (any(modp)) {		# multiple models - form table
	      opts <- dots[!modp]
	      mods <- c(list(object), dots[modp])
	      names(mods) <- sapply(as.list(mCall)[c(FALSE, TRUE, modp)],
				    as.character)
	      mods <- mods[order(sapply(lapply(mods, logLik, REML = FALSE),
					attr, "df"))]
	      calls <- lapply(mods, slot, "call")
	      data <- lapply(calls, "[[", "data")
	      if (any(data != data[[1]]))
		  stop("all models must be fit to the same data object")
	      header <- paste("Data:", data[[1]])
	      subset <- lapply(calls, "[[", "subset")
	      if (any(subset != subset[[1]]))
		  stop("all models must use the same subset")
	      if (!is.null(subset[[1]]))
		  header <-
		      c(header, paste("Subset", deparse(subset[[1]]), sep = ": "))
	      llks <- lapply(mods, logLik, REML = FALSE)
	      Df <- sapply(llks, attr, "df")
	      llk <- unlist(llks)
	      chisq <- 2 * pmax(0, c(NA, diff(llk)))
	      dfChisq <- c(NA, diff(Df))
	      val <- data.frame(Df = Df,
				AIC = sapply(llks, AIC),
				BIC = sapply(llks, BIC),
				logLik = llk,
				"Chisq" = chisq,
				"Chi Df" = dfChisq,
				"Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
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
            if (length(object@muEta))
              stop("single argument anova for GLMMs not yet implemented")
            if (length(object@V))
              stop("single argument anova for NLMMs not yet implemented")

            p <- object@dims[["p"]]
            ss <- (.Call(mer_update_projection, object)[[2]])^2
            names(ss) <- names(object@fixef)
            asgn <- attr(object@X, "assign")

            terms <- terms(object)
            nmeffects <- attr(terms, "term.labels")
            if ("(Intercept)" %in% names(ss))
              nmeffects <- c("(Intercept)", nmeffects)
            ss <- unlist(lapply(split(ss, asgn), sum))
            df <- unlist(lapply(split(asgn,  asgn), length))
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
      })

if (FALSE) {
    setMethod("confint", signature(object = "mer"),
              function(object, parm, level = 0.95, ...)
              .NotYetImplemented()
              )
}

setMethod("deviance", signature(object="mer"),
	  function(object, REML = NULL, ...)
      {
          if (missing(REML) || is.null(REML) || is.na(REML[1]))
	      REML <- object@dims[["REML"]]
	  object@deviance[[if(REML) "REML" else "ML"]]
      })

setMethod("fitted", signature(object = "mer"),
          function(object, ...)
          napredict(attr(object@frame, "na.action"), object@mu))

setMethod("formula", signature(x = "mer"),
	  function(x, ...)
	  x@call$formula
	  )

setMethod("logLik", signature(object="mer"),
	  function(object, REML = NULL, ...)
### Extract the log-likelihood or restricted log-likelihood
      {
          dims <- object@dims
          if (is.null(REML) || is.na(REML[1]))
              REML <- dims[["REML"]]
          val <- -deviance(object, REML = REML)/2
          attr(val, "nall") <- attr(val, "nobs") <- dims[["n"]]
          attr(val, "df") <-
              dims[["p"]] + dims[["np"]] + as.logical(dims[["useSc"]])
          attr(val, "REML") <-  as.logical(REML)
          class(val) <- "logLik"
          val
      })
setMethod("predict", signature(object="mer"),
          function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
                    interval = c("none", "confidence", "prediction"),
                    level = 0.95, type = c("response", "terms"),
                    terms = NULL, na.action = na.pass,
                    pred.var = res.var/weights, weights = 1, ...)
          {
            tt <- terms(object)
            if (missing(newdata) || is.null(newdata)) {
              mm <- X <- model.matrix(object)
              mmDone <- TRUE
              offset <- object@offset
            }
            else {
              Terms <- delete.response(tt)
              m <- model.frame(Terms, newdata, na.action = na.action,
                               xlev = object$xlevels)
              if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, m)
              X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
              offset <- if (!is.null(off.num <- attr(tt, "offset")))
                eval(attr(tt, "variables")[[off.num + 1]], newdata)
              else if (!is.null(object$offset))
                eval(object$call$offset, newdata)
              mmDone <- FALSE
            }
            n <- length(residuals(object))
            predictor <- drop(X %*% fixef(object))
            if (length(offset))
              predictor <- predictor + offset
            return(predictor)

            interval <- match.arg(interval)
            if (interval == "prediction") {
                if (missing(newdata))
                    warning("Predictions on current data refer to _future_ responses\n")
                if (missing(newdata) && missing(weights)) {
                    w <- weights.default(object)
                    if (!is.null(w)) {
                        weights <- w
                        warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
                    }
                }
                if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
                    missing(pred.var))
                    warning("Assuming constant prediction variance even though model fit is weighted\n")
                if (inherits(weights, "formula")) {
                    if (length(weights) != 2L)
                        stop("'weights' as formula should be one-sided")
                    d <- if (missing(newdata) || is.null(newdata))
                        model.frame(object)
                    else newdata
                    weights <- eval(weights[[2L]], d, environment(weights))
                }
            }
            type <- match.arg(type)
            if (se.fit || interval != "none") {
                res.var <- if (is.null(scale)) {
                    r <- object$residuals
                    w <- object$weights
                    rss <- sum(if (is.null(w)) r^2 else r^2 * w)
                    df <- n - p
                    rss/df
                }
                else scale^2
                if (type != "terms") {
                    if (p > 0) {
                        XRinv <- if (missing(newdata) && is.null(w))
                            qr.Q(object$qr)[, p1, drop = FALSE]
                        else X[, piv] %*% qr.solve(qr.R(object$qr)[p1,
                                                                   p1])
                        ip <- drop(XRinv^2 %*% rep(res.var, p))
                    }
                    else ip <- rep(0, n)
                }
            }
            if (type == "terms") {
                if (!mmDone) {
                    mm <- model.matrix(object)
                    mmDone <- TRUE
                }
                aa <- attr(mm, "assign")
                ll <- attr(tt, "term.labels")
                hasintercept <- attr(tt, "intercept") > 0L
                if (hasintercept)
                    ll <- c("(Intercept)", ll)
                aaa <- factor(aa, labels = ll)
                asgn <- split(order(aa), aaa)
                if (hasintercept) {
                    asgn$"(Intercept)" <- NULL
                    if (!mmDone) {
                        mm <- model.matrix(object)
                        mmDone <- TRUE
                    }
                    avx <- colMeans(mm)
                    termsconst <- sum(avx[piv] * beta[piv])
                }
                nterms <- length(asgn)
                if (nterms > 0) {
                    predictor <- matrix(ncol = nterms, nrow = NROW(X))
                    dimnames(predictor) <- list(rownames(X), names(asgn))
                    if (se.fit || interval != "none") {
                        ip <- matrix(ncol = nterms, nrow = NROW(X))
                        dimnames(ip) <- list(rownames(X), names(asgn))
                        Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
                    }
                    if (hasintercept)
                        X <- sweep(X, 2L, avx, check.margin = FALSE)
                    unpiv <- rep.int(0L, NCOL(X))
                    unpiv[piv] <- p1
                    for (i in seq.int(1L, nterms, length.out = nterms)) {
                        iipiv <- asgn[[i]]
                        ii <- unpiv[iipiv]
                        iipiv[ii == 0L] <- 0L
                        predictor[, i] <- if (any(iipiv > 0L))
                            X[, iipiv, drop = FALSE] %*% beta[iipiv]
                        else 0
                        if (se.fit || interval != "none")
                            ip[, i] <- if (any(iipiv > 0L))
                                as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,
                                                     , drop = FALSE])^2 %*% rep.int(res.var,
                                                       p)
                            else 0
                    }
                    if (!is.null(terms)) {
                        predictor <- predictor[, terms, drop = FALSE]
                        if (se.fit)
                            ip <- ip[, terms, drop = FALSE]
                    }
                }
                else {
                    predictor <- ip <- matrix(0, n, 0)
                }
                attr(predictor, "constant") <- if (hasintercept)
                    termsconst
                else 0
            }
            if (interval != "none") {
                tfrac <- qt((1 - level)/2, df)
                hwid <- tfrac * switch(interval, confidence = sqrt(ip),
                                       prediction = sqrt(ip + pred.var))
                if (type != "terms") {
                    predictor <- cbind(predictor, predictor + hwid %o%
                                       c(1, -1))
                    colnames(predictor) <- c("fit", "lwr", "upr")
                }
                else {
                    lwr <- predictor + hwid
                    upr <- predictor - hwid
                }
            }
            if (se.fit || interval != "none")
                se <- sqrt(ip)
            if (missing(newdata) && !is.null(na.act <- object$na.action)) {
                predictor <- napredict(na.act, predictor)
                if (se.fit)
                    se <- napredict(na.act, se)
            }
            if (type == "terms" && interval != "none") {
                if (missing(newdata) && !is.null(na.act)) {
                    lwr <- napredict(na.act, lwr)
                    upr <- napredict(na.act, upr)
                }
                list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
                     df = df, residual.scale = sqrt(res.var))
            }
            else if (se.fit)
                list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
            else predictor
        })



setMethod("residuals", signature(object = "mer"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

setMethod("resid", signature(object = "mer"),
	  function(object, ...)
          napredict(attr(object@frame, "na.action"), object@resid))

setMethod("simulate", "mer",
          function(object, nsim = 1, seed = NULL, ...)
      {
	  if(!is.null(seed)) set.seed(seed)
	  if(!exists(".Random.seed", envir = .GlobalEnv))
	      runif(1)		     # initialize the RNG if necessary
          RNGstate <- .Random.seed
          ## FIXME: implement offset, should be fairly easy (?)
          dims <- object@dims
          sigma <- sigma(object)
          etasim.fix <- as.vector(object@X %*% fixef(object))   # fixed-effect contribution
          if (length(offset <- object@offset)>0) {
            etasim.fix <- etasim.fix+offset
          } 
          etasim.reff <- as(t(object@A) %*%    # UNSCALED random-effects contribution
                            matrix(rnorm(nsim * dims[["q"]]), nc = nsim),
                            "matrix")
          if (length(object@V) == 0 && length(object@muEta) == 0) {
            etasim.resid <- matrix(rnorm(nsim * dims["n"]), nc = nsim) ## UNSCALED residual
            etasim <- etasim.fix + sigma*(etasim.reff+etasim.resid)
            return(drop(etasim))
          }
          if (length(object@muEta)>0) {
      ## GLMM
      ## n.b. DON'T scale random-effects
      ## if sigma!=1, it applies to the "quasi"- part of the model
      etasim <- etasim.fix+etasim.reff 
      family <- object@call$family
      if(is.symbol(family)) family <- as.character(family)
      if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
      if(is.function(family)) family <- family()
      if(is.null(family$family)) stop("'family' not recognized")
      musim <- family$linkinv(etasim)
      n <- length(musim) ## FIXME: or could be dims["n"]?
      vsim <- switch(family$family,
                     poisson=rpois(n,lambda=musim),
                     binomial={
                       resp <- model.response(object@frame)
                       bernoulli <- !is.matrix(resp)
                       if (bernoulli) {
                         rbinom(n,prob=musim,size=1)
                       } else {
                         nresp <- nrow(resp)
                         ## FIXME: should "N-size" (column 2) be named?
                         ## copying structures from stats/R/family.R
                         sizes <- rowSums(resp)
                         Y <- rbinom(n, size = sizes, prob = musim)
                         YY <- cbind(Y, sizes - Y)
                         yy <- lapply(split(YY,gl(nsim,nresp,2*nsim*nresp)),
                                matrix,ncol=2,dimnames=list(NULL,colnames(resp)))
                         ## colnames() <- colnames(resp)
                         ## yy <- split(as.data.frame(YY),
                         ## rep(1:nsim,each=length(sizes)))
                         names(yy) <- paste("sim",seq_along(yy),sep="_")
                         yy
                         ## yy <- as.data.frame(yy)
                       }
                     },
                     stop("simulation not implemented for family",
                          family$family))
                       
    }
    if (!(family$family=="binomial" && !bernoulli)) {
      vsim <- matrix(vsim,nc=nsim)
      return(drop(vsim))
    } else if (nsim==1) return(vsim[[1]]) else return(vsim)
    stop("simulate method for NLMMs not yet implemented")
})

setMethod("summary", signature(object = "mer"),
	  function(object, ...)
      {
          REML <- object@dims[["REML"]]
          fcoef <- fixef(object)
          vcov <- vcov(object)
          corF <- vcov@factors$correlation
          dims <- object@dims
          coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
          llik <- logLik(object, REML)
          dev <- object@deviance
          mType <- if((non <- as.logical(length(object@V)))) "NMM" else "LMM"
          if (gen <- as.logical(length(object@muEta)))
              mType <- paste("G", mType, sep = '')
          mName <- switch(mType, LMM = "Linear", NMM = "Nonlinear",
                          GLMM = "Generalized linear",
                          GNMM = "Generalized nonlinear")
	  method <- {
	      if (mType == "LMM")
		  if(REML) "REML" else "maximum likelihood"
	      else
		  paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
			"adaptive Gaussian Hermite",
			"approximation")
	  }
          AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                 logLik = as.vector(llik),
                                 deviance = dev[["ML"]],
                                 REMLdev = dev[["REML"]],
                                 row.names = "")
          if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
          varcor <- VarCorr(object)
          REmat <- formatVC(varcor)
          if (is.na(attr(varcor, "sc")))
              REmat <- REmat[-nrow(REmat), , drop = FALSE]

          if (nrow(coefs) > 0) {
              if (!dims[["useSc"]]) {
                  coefs <- coefs[, 1:2, drop = FALSE]
                  stat <- coefs[,1]/coefs[,2]
                  pval <- 2*pnorm(abs(stat), lower = FALSE)
                  coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
              } else {
                  stat <- coefs[,1]/coefs[,2]
                  ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
                  coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
              }
          } ## else : append columns to 0-row matrix ...
          new("summary.mer",
              object,
              methTitle = paste(mName, "mixed model fit by", method),
              logLik = llik,
              ngrps = sapply(object@flist, function(x) length(levels(x))),
              sigma = sigma(object),
              coefs = coefs,
              vcov = vcov,
              REmat = REmat,
              AICtab= AICframe
              )
      })## summary()

setMethod("model.frame", signature(formula = "mer"),
	  function(formula, ...) formula@frame)

setMethod("model.matrix", signature(object = "mer"),
	  function(object, ...) object@X)

setMethod("terms", signature(x = "mer"),
	  function(x, ...) attr(x@frame, "terms"))

setMethod("update", signature(object = "mer"),
	  function(object, formula., ..., evaluate = TRUE)
      {
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
      })

setMethod("vcov", signature(object = "mer"),
	  function(object, ...)
### Extract the conditional variance-covariance matrix of the fixed effects
      {
          rr <- as(sigma(object)^2 *
                   chol2inv(object@RX, size = object@dims['p']), "dpoMatrix")
          nms <- colnames(object@X)
          dimnames(rr) <- list(nms, nms)
          rr@factors$correlation <- as(rr, "corMatrix")
          rr
      })

setMethod("with", signature(data = "mer"),
	  function(data, expr, ...) {
	      dat <- eval(data@call$data)
	      if (!is.null(na.act <- attr(data@frame, "na.action")))
		  dat <- dat[-na.act, ]
	      lst <- c(list(. = data), data@flist, data@frame, dat)
	      eval(substitute(expr), lst[unique(names(lst))])
	  })

### Show and print methods and utilities for them

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x, maxlen) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }, maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
	cbind(reMat, rBind(corr, rep.int("", ncol(corr))))
    } else reMat
}

## This is modeled a bit after  print.summary.lm :
printMer <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = TRUE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    REML <- so@dims[["REML"]]
    llik <- so@logLik
    dev <- so@deviance
    dims <- x@dims

    cat(so@methTitle, "\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:", deparse(x@call$subset),"\n")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    if (is.na(so@sigma))
	cat("\nEstimated scale (compare to 1):",
            sqrt(exp(so@deviance[["lr2"]])/so@dims[["n"]]), "\n")
    if (nrow(so@coefs) > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(correlation) {
	    corF <- so@vcov@factors$correlation
	    if (!is.null(corF)) {
		p <- ncol(corF)
		if (p > 1) {
		    rn <- rownames(so@coefs)
		    rns <- abbreviate(rn, minlen=11)
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			corf <- as(corF, "matrix")
			dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			print(symnum(corf))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p,
                                       dimnames = list(rns, abbreviate(rn, minlen=6)))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}

setMethod("print", "mer", printMer)
setMethod("show", "mer", function(object) printMer(object))

printNlmer <- function(x, digits = max(3, getOption("digits") - 3),
                       correlation = TRUE, symbolic.cor = FALSE,
                       signif.stars = getOption("show.signif.stars"), ...)
### FIXME: Does nlmer need a separate show method?
{
    dims <- x@dims
    cat("Nonlinear mixed model fit by Laplace\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:", deparse(x@call$subset),"\n")

    cat("Random effects:\n")
    print(formatVC(VarCorr(x)), quote = FALSE,
          digits = max(3, getOption("digits") - 3))

    cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
    ngrps <- sapply(x@flist, function(x) length(levels(x)))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    cat("\nFixed effects:\n")
    print(x@fixef)
    invisible(x)
}

setMethod("refit", signature(object = "mer", newresp = "numeric"),
          function(object, newresp, ...)
      {
          newresp <- as.double(newresp[!is.na(newresp)])
          stopifnot(length(newresp) == object@dims[["n"]])
          object@y <- newresp
          mer_finalize(object)
      })

## Contributed by Ben Bolker
setMethod("refit", signature(object = "mer", newresp = "matrix"),
          function(object, newresp, ...)
      {
          stopifnot(ncol(newresp) == 2,
                    all(!is.na(wts <- rowSums(newresp))),
                    length(wts) == object@dims[["n"]])
          object@y <- newresp[,1]/wts
          object@pWt <- wts
          mer_finalize(object)
      })


BlockDiagonal <- function(lst)
{
    stopifnot(is(lst, "list"))
    lst <- lapply(lapply(lst, as, Class = "generalMatrix"),
                  as, Class = "TsparseMatrix")
    isSquare <- function(x) nrow(x) == ncol(x)
    stopifnot(all(sapply(lst, isSquare)),
              all(sapply(lst, is, class2 = "dMatrix")))
    if ((nl <- length(lst)) == 1) return(lst[[1]])

    offsets <- c(0L, cumsum(sapply(lst, ncol)))
    new("dgTMatrix", Dim = rep.int(offsets[nl + 1], 2),
        i = unlist(lapply(1:nl, function(i) lst[[i]]@i + offsets[i])),
        j = unlist(lapply(1:nl, function(i) lst[[i]]@j + offsets[i])),
        x = unlist(lapply(lst, slot, "x")))
}

setMethod("expand", signature(x = "mer"),
          function(x, sparse = TRUE, ...)
      {
          ind <- seq_along(ST <- x@ST)

          if (!sparse) {
              elexpand <- function(mat)
                  list(T = new("dtrMatrix", uplo = "L", diag = "U",
                       x = as.vector(mat),
                       Dim = dim(mat), Dimnames = dimnames(mat)),
                       S = Diagonal(x = diag(mat)))
              fl <- x@flist
              if (all(attr(fl, "assign") == ind))
                  names(ST) <- names(fl)
              return(lapply(ST, elexpand))
          }

          ## 'Sparse' case :

          nc <- sapply(ST, ncol)
          if(!all(nc >= 1)) stop("some ST entries lack  ncol(.) >= 1")
          nlev <- diff(x@Gp) %/% nc

          Sblock <- function(i) rep(diag(ST[[i]]), each = nlev[i])
          Smat <- Diagonal(x = unname(unlist(lapply(ind, Sblock))))
          if (max(nc) == 1) {
              Tmat <- Matrix:::.diag2tT(Diagonal(ncol(Smat)), uplo="L")
          } else {
              Tblock <- function(i)
              {
                  if (nc[i] == 1) return(Matrix:::.diag2tT(Diagonal(nlev[i]), uplo="L"))
                  STi <- ST[[i]]
                  nci <- nc[i]
                  lt <- lower.tri(STi)
                  offsets <- (1:nlev[i]) - 1L
                  ij <- nlev[i] * (which(lt, arr.ind=TRUE) - 1)
                  new("dtTMatrix", Dim = rep.int(nlev[i] * nci, 2), uplo = "L", diag = "U",
                      x = rep(STi[lt], each = nlev[i]),
                      i = as.integer(outer(offsets, ij[, 1], "+")),
                      j = as.integer(outer(offsets, ij[, 2], "+")))
              }
              Tmat <- as(as(BlockDiagonal(lapply(ind, Tblock)), "triangularMatrix"),
                         "CsparseMatrix")
          }
          list(sigma = sigma(x), P = as(x@L@perm + 1L, "pMatrix"),
               T = as(Tmat, "CsparseMatrix"), S = Smat)
      })

#### Methods for secondary, derived classes

setMethod("deviance", signature(object = "summary.mer"), function(object) object@deviance)
setMethod("logLik", signature(object = "summary.mer"), function(object) object@logLik)
setMethod("vcov", signature(object = "summary.mer"), function(object) object@vcov)
setMethod("summary", signature(object = "summary.mer"), function(object) object)

#### Methods to produce specific plots

plot.coef.mer <- function(x, y, ...)
{
    varying <- unique(do.call("c",
                              lapply(x, function(el)
                                     names(el)[sapply(el,
                                                      function(col)
                                                      any(col != col[1]))])))
    gf <- do.call("rBind", lapply(x, "[", j = varying))
    gf$.grp <- factor(rep(names(x), sapply(x, nrow)))
    switch(min(length(varying), 3),
           qqmath(eval(substitute(~ x | .grp,
                                  list(x = as.name(varying[1])))), gf, ...),
           xyplot(eval(substitute(y ~ x | .grp,
                                  list(y = as.name(varying[1]),
                                       x = as.name(varying[2])))), gf, ...),
           splom(~ gf | .grp, ...))
}

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

qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        hw <- 1.96 * se
        list(ylim = range(y - hw, y + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
        panel.grid(h = -1,v = -1)
        panel.abline(h = 0)
        x <- as.numeric(x)
        y <- as.numeric(y)
        se <- as.numeric(se[subscripts])
        ly <- y - 1.96 * se
        uy <- y + 1.96 * se
        panel.segments(x, y - 1.96*se, x, y + 1.96 * se,
                       col = 'black')
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
            xyplot(unlist(x)[ord] ~
                   rep(qnorm((rr - 0.5)/nr), ncol(x)) | ind[ord],
                   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, aspect = 1, ...)
        } else {
            qqmath(~values|ind, stack(x),
                   scales = list(y = list(relation = "free")),
                   xlab = "Standard normal quantiles",
                   ylab = NULL, ...)
        }
    }
    lapply(x, f)
}

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

#### Creating and displaying a Markov Chain Monte Carlo sample from
#### the posterior distribution of the parameters

setMethod("mcmcsamp", signature(object = "mer"),
	  function(object, n = 1, verbose = FALSE, saveb = FALSE, ...)
### Generate a Markov chain Monte Carlo sample from the posterior distribution
### of the parameters in a linear mixed model
      {
          object@fixef <- fixef(object) # force a copy
          n <- max(1, as.integer(n)[1])
          dd <- object@dims
          ranef <- matrix(numeric(0), nrow = dd[["q"]], ncol = 0)
          if (saveb) ranef <- matrix(object@ranef, nrow = dd[["q"]], ncol = n)
          sigma <- matrix(unname(sigma(object)), nrow = 1,
                          ncol = (if (dd[["useSc"]]) n else 0))
          ff <- object@fixef
          fixef <- matrix(ff, dd[["p"]], n)
          rownames(fixef) <- names(ff)
          ans <- new("merMCMC",
                     Gp = object@Gp,
                     ST = matrix(.Call(mer_ST_getPars, object), dd[["np"]], n),
                     call = object@call,
                     dims = object@dims,
                     deviance = rep(unname(object@deviance[["ML"]]), n),
                     fixef = fixef,
                     nc = sapply(object@ST, nrow),
                     ranef = ranef,
                     sigma = sigma)
          .Call(mer_MCMCsamp, ans, object)
      })

setMethod("HPDinterval", signature(object = "merMCMC"),
          function(object, prob = 0.95, ...)
      {
          nms <- c("fixef", "ST")
          if (length(object@sigma)) nms <- c(nms, "sigma")
          if (length(object@ranef)) nms <- c(nms, "ranef")
          names(nms) <- nms
          lapply(lapply(nms, slot, object = object),
                 HPDinterval, prob = prob)
      })

setMethod("HPDinterval", signature(object = "matrix"),
          function(object, prob = 0.95, ...)
      {
          if (ncol(object) > nrow(object))
              object <- t(object)
          vals <- apply(object, 2, sort)
          if (!is.matrix(vals))
              stop("object must have nsamp > 1")
          nsamp <- nrow(vals)
          npar <- ncol(vals)
          gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
          init <- 1:(nsamp - gap)
          inds <- apply(vals[init + gap, , drop = FALSE] -
                        vals[init, , drop = FALSE], 2, which.min)
          ans <- cbind(vals[cbind(inds, 1:npar)],
                       vals[cbind(inds + gap, 1:npar)])
          dimnames(ans) <- list(colnames(object), c("lower", "upper"))
          attr(ans, "Probability") <- gap/nsamp
          ans
      })

### FIXME: Watch the names of the variance components here
setMethod("VarCorr", signature(x = "merMCMC"),
          function(x, type = c("raw", "varcov", "sdcorr", "logs"), ...)
      {
          if ("raw" == (type <- match.arg(type))) {
              ST <- t(x@ST)
              colnames(ST) <- paste("ST", 1:ncol(ST), sep = '')
              if (length(x@sigma)) return(cbind(ST, sigma = as.vector(x@sigma)))
              return(ST)
          }
          .Call(merMCMC_VarCorr, x, match(type, c("raw", "varcov", "sdcorr", "logs")))
      })

setMethod("as.matrix", signature(x = "merMCMC"),
          function(x, ...)
          cbind(t(x@fixef), VarCorr(x, ...)))

setMethod("as.data.frame", signature(x = "merMCMC"),
          function(x, row.names = NULL, optional = FALSE, ...)
          as.data.frame(as.matrix(x, ...), row.names = row.names, optional = optional, ...))

setAs("merMCMC", "data.frame", function(from) as.data.frame(from))

aslatticeframe <- function(x, ...)
{
    fr <- as.data.frame(x, ...)
    data.frame(dat = unlist(fr),
               par = gl(ncol(fr), nrow(fr), labels = colnames(fr)),
               iter = rep(1:nrow(fr), ncol(fr)))
}

## FIXME: More care should be taken to avoid duplicate argument names
## in the eventual call to lattice functions. Accumulate the arguments
## in a list and use do.call instead of direct calls.

setMethod("xyplot", signature(x = "merMCMC"),
          function(x, data, ...)
      {
          pfr <- aslatticeframe(x, ...)
          xyplot(dat ~ iter|par, pfr,
                 xlab = "Iteration number", ylab = NULL,
                 scales = list(x = list(axs = 'i'),
                 y = list(relation = "free", rot = 0)),
                 type = c("g", "l"),
                 layout = c(1, length(levels(pfr$par))),
                 strip = FALSE, strip.left = TRUE, ...)
      })

setMethod("densityplot", signature(x = "merMCMC"),
          function(x, data, ...)
          densityplot(~ dat | par, aslatticeframe(x, ...),
                      scales = list(relation = 'free'), ...)
          )

setMethod("qqmath", signature(x = "merMCMC"),
          function(x, data, ...)
          qqmath(~ dat | par, aslatticeframe(x, ...),
                 scales = list(y = list(relation = 'free')), ...)
          )


abbrvNms <- function(gnm, cnms)
### Abbreviate names of columns in grouping factors
### gnm - group name
### cnms - column names
{
    ans <- paste(abbreviate(gnm), abbreviate(cnms), sep = '.')
    if (length(cnms) > 1) {
	anms <- lapply(cnms, abbreviate, minlength = 3)
	nmmat <- outer(anms, anms, paste, sep = '.')
	ans <- c(ans, paste(abbreviate(gnm, minlength = 3),
			    nmmat[upper.tri(nmmat)], sep = '.'))
    }
    ans
}

mcmccompnames <- function(ans, object, saveb, trans, glmer, deviance)
### Mangle the names of the columns of the mcmcsamp result ans
### This operation is common to the methods for "lmer" and "glmer"
{
    gnms <- names(object@flist)
    cnms <- lapply(object@ST, colnames)
    ff <- fixef(object)
    colnms <- c(names(ff), if (glmer) character(0) else "sigma^2",
                unlist(lapply(seq_along(gnms),
                              function(i)
                              abbrvNms(gnms[i],cnms[[i]]))))
    if (trans) {
        ## parameter type: 0 => fixed effect, 1 => variance,
        ##		 2 => covariance
        ptyp <- c(integer(length(ff)), if (glmer) integer(0) else 1:1,
                  unlist(lapply(seq_along(gnms),
                                function(i)
                            {
                                k <- length(cnms[[i]])
                                rep(1:2, c(k, (k*(k-1))/2))
                            })))
        colnms[ptyp == 1] <-
            paste("log(", colnms[ptyp == 1], ")", sep = "")
        colnms[ptyp == 2] <-
            paste("atanh(", colnms[ptyp == 2], ")", sep = "")
    }
    if (deviance) colnms <- c(colnms, "deviance")
### FIXME: this will fail for a mer2 object
    if(saveb) {## maybe better colnames, "RE.1","RE.2", ... ?
        .NotYetImplemented()
        rZy <- object@rZy
        colnms <- c(colnms,
                    paste("b", sprintf(paste("%0",
                                             1+floor(log(length(rZy),10)),
                                             "d", sep = ''),
                                       seq_along(rZy)),
                          sep = '.'))
    }
    colnames(ans) <- colnms
    ans
}

devvals <- function(fm, pmat, sigma1 = FALSE)
{
    if (!is(fm, "mer"))
        stop('fm must be an "mer" fitted model')
### FIXME: add a check in here for glmer and nlmer
    np <- length(p0 <- .Call(mer_ST_getPars, fm))
    pmat <- as.matrix(pmat)
    if (ncol(pmat) != np + sigma1)
        stop(gettextf("pmat must have %d columns", np + sigma1))
    storage.mode(pmat) <- "double"
    if (is.null(pnms <- dimnames(pmat)[[2]]))
        pnms <- c(if(sigma1) character(0) else "sigma",
                  paste("th", seq_len(np), sep = ""))
    dev <- fm@deviance
    ans <- matrix(0, nrow(pmat), ncol(pmat) + length(dev),
                  dimnames = list(NULL, c(pnms, names(dev))))
    for (i in seq_len(nrow(pmat))) {
        .Call(mer_ST_setPars, fm,
              ## This expression does not allow for correlated random
              ## effects.  It would be best to make the appropriate
              ## changes in the C code for ST_setPars.
              if (sigma1) pmat[i,-1]/pmat[i,1] else pmat[i,])
        ans[i, ] <- c(pmat[i, ], .Call(mer_update_RX, fm))
    }
    .Call(mer_ST_setPars, fm, p0)
    .Call(mer_update_RX, fm)
    as.data.frame(ans)
}

#### Odds and ends

## simulestimate <- function(x, FUN, nsim = 1, seed = NULL, control = list())
## {
##     FUN <- match.fun(FUN)
##     stopifnot((nsim <- as.integer(nsim[1])) > 0,
## 	      inherits(x, "lmer"))
##     if (!is.null(seed)) set.seed(seed)
##     ## simulate the linear predictors
##     lpred <- .Call(mer_simulate, x, nsim)
##     sc <- abs(x@devComp[8])
##     ## add fixed-effects contribution and per-observation noise term
##     lpred <- lpred + drop(x@X %*% fixef(x)) + rnorm(prod(dim(lpred)), sd = sc)

##     cv <- do.call(lmerControl, control)
##     Omega <- x@Omega
##     x@wrkres <- x@y <- lpred[,1]
##     .Call(mer_update_ZXy, x)
##     LMEoptimize(x) <- cv
##     template <- FUN(x)
##     if (!is.numeric(template))
##         stop("simulestimate currently only handles functions that return numeric vectors")
##     ans <- matrix(template, nr = nsim, nc = length(template), byrow = TRUE)
##     colnames(ans) <- names(template)
##     for (i in 1:nsim) {
##         x@wrkres <- x@y <- lpred[,i]
##         x@Omega <- Omega
##         .Call(mer_update_ZXy, x)
##         LMEoptimize(x) <- cv
##         foo <- try(FUN(x))
##         ans[i,] <- if (inherits(foo, "try-error")) NA else foo
##     }
##     ans
## }

hatTrace <- function(x)
{
    .NotYetImplemented()
    stopifnot(is(x, "mer"))
}

ST2Omega <- function(ST)
### Temporary function to convert the ST representation of the
### relative variance-covariance matrix returned by lmer into the
### Omega representation required by lmer
{
    if (nrow(ST) == 1) return(as(1/ST^2, "dpoMatrix"))
    dd <- diag(ST)
    T <- as(ST, "dtrMatrix")
    T@diag <- "U"
    crossprod(solve(T)/dd)
}


## setMethod("simulate", signature(object = "mer"),
## 	  function(object, nsim = 1, seed = NULL, ...)
##       {
## 	  if(!exists(".Random.seed", envir = .GlobalEnv))
## 	      runif(1)		     # initialize the RNG if necessary
## 	  if(is.null(seed))
## 	      RNGstate <- .Random.seed
## 	  else {
## 	      R.seed <- .Random.seed
## 	      set.seed(seed)
## 	      RNGstate <- structure(seed, kind = as.list(RNGkind()))
## 	      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
## 	  }

##           stopifnot((nsim <- as.integer(nsim[1])) > 0,
##                     inherits(object, "lmer"))
## 	  ## similate the linear predictors
## 	  lpred <- .Call(mer_simulate, object, nsim)
## 	  sc <- abs(object@devComp[8])

## 	  ## add fixed-effects contribution and per-observation noise term
## 	  lpred <- as.data.frame(lpred + drop(object@X %*% fixef(object)) +
## 				 rnorm(prod(dim(lpred)), sd = sc))
## 	  ## save the seed
## 	  attr(lpred, "seed") <- RNGstate
## 	  lpred
##       })

## We need to define an S4 print method, since using an S3 print
## method fails as soon as you call print() explicitly, e.g. when
## wanting to specify options.

## calculates degrees of freedom for fixed effects Wald tests
## This is a placeholder.  The answers are generally wrong.  It will
## be very tricky to decide what a 'right' answer should be with
## crossed random effects.

## setMethod("getFixDF", signature(object="mer"),
## 	  function(object, ...) {
## 	      devc <- object@devComp
## 	      rep(as.integer(devc[1]- devc[2]), devc[2])
## 	  })

## simss <- function(fm0, fma, nsim)
## {
##     ysim <- simulate(fm0, nsim)
##     cv <- list(gradient = FALSE, msMaxIter = 200:200,
## 	       msVerbose = 0:0)
##     sapply(ysim, function(yy) {
## 	.Call(mer_update_y, fm0, yy)
## 	LMEoptimize(fm0) <- cv
## 	.Call(mer_update_y, fma, yy)
## 	LMEoptimize(fma) <- cv
## 	exp(c(H0 = fm0@devComp[[["logryy2"]]],
## 	      Ha = fma@devComp[[["logryy2"]]]))
##     })
## }

## setMethod("denomDF", "mer",
##           function(x, ...)
##       {
##           mm <- x@X
##           aa <- attr(mm, "assign")
##           tt <- x@terms
##           if (!isNested(x))
##               return(list(coef = as.numeric(rep(NA, length(x@fixef))),
##                           terms = as.numeric(rep(NA,
##                           length(attr(tt, "order"))))))
##           hasintercept <- attr(tt, "intercept") > 0
##           ## check which variables vary within levels of grouping factors
##           vars <- eval(attr(tt, "variables"), x@frame)
##           fl <- x@flist
##           vv <- matrix(0:0, nrow = length(vars), ncol = length(fl),
##                         dimnames = list(NULL, names(fl)))
##           ## replace this loop by C code.
##           for (i in 1:nrow(ans))        # check if variables vary within factors
##               for (j in 1:ncol(ans))
##                   ans[i,j] <- all(tapply(vars[[i]], fl[[j]],
##                                          function(x) length(unique(x)) == 1))
##           ## which terms vary within levels of which grouping factors?
##           tv <- crossprod(attr(tt, "factors"), !ans)
##           ## maximum level at which the term is constant
##           ml <- apply(tv, 1, function(rr) max(0, which(as.logical(rr))))
##           ## unravel assignment applied to terms
##           ll <- attr(tt, "term.labels")
##           if (hasintercept)
##               ll <- c("(Intercept)", ll)
##           aaa <- factor(aa, labels = ll)
##           asgn <- split(order(aa), aaa)
##           nco <- lapply(asgn, length)   # number of coefficients per term
##           nlev <- lapply(fl, function(x) length(levels(x)))
##           if (hasintercept) asgn$"(Intercept)" <- NULL
##           list(ml = ml, nco = nco, nlev = nlev)
##       })

## Utilities for the fitted mer object
slotsz <- function(obj)
    rev(sort(sapply(slotNames(obj), function(s) object.size(slot(obj, s)))))

slotApply <- function(object, f, ..., simplify = FALSE) {
   .localFun <- function(what, ...) f(slot(object, what), ...)
   sapply(slotNames(object), .localFun, ..., simplify = simplify)
}


yfrm <- function(fm)
{
    stopifnot(is(fm, "mer"))
    snr <- slotApply(fm, function(x)
                 {
                     if (is(x, "matrix") ||
                         is(x, "data.frame") ||
                         is(x, "numeric")) return (NROW(x))
                     0
                 }, simplify = TRUE)
    snr <- snr[snr > 0 & !(names(snr) %in%
                           c("Gp", "dims", "deviance", "frame", "flist", "X"))]
    fr <- cbind(fm@frame, fm@flist[1:NROW(fm@frame), !(names(fm@flist) %in%
                                     names(fm@frame))])
    n <- NROW(fr)
    if (NROW(fm@X) == n)
        fr <- cbind(fr, X = fm@X, Xbeta = fm@X %*% fm@fixef,
                    Zb = crossprod(fm@Zt, fm@ranef)@x)
    do.call(cbind, c(list(fr), sapply(names(which(snr == NROW(fr))),
                                      slot, object = fm, simplify = FALSE)))
}

##' Evaluate conditional components of an LMM.
##'
##' Evaluate conditional components of a linear mixed model for a grid of ST
##' parameter values.
##'
##' @param fm - a fitted linear mixed model
##' @param parmat - a numeric matrix whose rows constitute suitable parameter
##'     values for fm@ST
##' @param type - which slot to extract
##' @return a data frame of deviance values or fixed-effects or random effects
##' @keywords models
devmat <-
    function(fm, parmat, slotname = c("deviance", "fixef", "ranef", "u"), ...)
{
    stopifnot(is(fm, "mer"))
    dd <- fm@dims
    stopifnot(dd[["fTyp"]] == 2L, # gaussian family
              dd[["lTyp"]] == 5L, # identity link
              dd[["vTyp"]] == 1L, # variance function is "constant"
              length(fm@V) == 0L, # nonlinear parameter gradient is identity
              length(fm@muEta) == 0L) # eta -> mu map is identity
    oldpars <- .Call(mer_ST_getPars, fm)

    parmat <- as.matrix(parmat)
    storage.mode(parmat) <- "double"
    if (ncol(parmat) == dd[["np"]])
        parmat <- t(parmat)             # parameter vectors as columns
    stopifnot(nrow(parmat) == dd[["np"]])
    slotname <- match.arg(slotname)

    slotval <- function(x) {            # function to apply
        .Call(mer_ST_setPars, fm, x)
        .Call(mer_update_dev, fm)
        if (slotname != "deviance") .Call(mer_update_ranef, fm)
        slot(fm, slotname)
    }
    ans <- apply(parmat, 2, slotval)
    slotval(oldpars)                    # restore the fitted model
    as.data.frame(t(rbind(parmat, ans)))
}

##' Find terms associated with grouping factor names.

##' Determine the random-effects associated with particular grouping
##' factors.

##' @param fm a fitted model object of S4 class "mer"
##' @param fnm one or more grouping factor names, as a character vector

##' @return a list of indices of terms
##' @keywords models
##' @export
##' @examples
##' fm1 <- lmer(strength ~ (1|batch) + (1|sample), Pastes)
##' whichterms(fm1)
whichterms <- function(fm, fnm = names(fm@flist))
{
    stopifnot(is(fm, "mer"), is.character(fnm))
    fl <- fm@flist
    asgn <- attr(fl, "assign")
    fnms <- names(fl)
    stopifnot(all(fnm %in% fnms))
    if (is.null(names(fnm))) names(fnm) <- fnm

    lapply(fnm, function(nm) which(asgn == match(nm, fnms)))
}

##' Random-effects indices by term

##' Returns a list of indices into the ranef vector by random-effects
##' terms.

##' @param Gp the Gp slot from an mer object

##' @return a list of random-effects indices
##' @keywords models
reinds <- function(Gp)
{
    lens <- diff(Gp)
    lapply(seq_along(lens), function(i) Gp[i] + seq_len(lens[i]))
}

##' Random-effects indices associated with grouping factor names

##' Determine the random-effects indices with particular grouping
##' factors.

##' @param fm a fitted model object of S4 class "mer"
##' @param fnm one or more grouping factor names, as a character vector

##' @return a list of indices of terms
##' @keywords models
##' @export
##' @examples
##' fm1 <- lmer(strength ~ (1|batch) + (1|sample), Pastes)
##' whichreind(fm1)
whichreind <- function(fm, fnm = names(fm@flist))
    lapply(whichterms(fm, fnm),
           function (ind) unlist(reinds(fm@Gp)[ind]))


## For Matrix API change (Oct.2009) - silence the warning:
## don't -- we do not call it in this version of lme4:
## assign("det_CHMfactor.warn", FALSE, envir = Matrix:::.MatrixEnv)
