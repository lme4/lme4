## Multilevel Item Response Theory model

mlirt <-
    function(formula, data, family = NULL, method = c("REML", "ML"),
             control = list(), start = NULL, verbose = FALSE,
             subset, weights, na.action, offset, contrasts = NULL,
             model = TRUE, x = TRUE, ...)
{
    mc <- match.call()
                                        # Evaluate and check the family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    stopifnot(length(formula <- as.formula(formula)) == 3,
              all.equal(family$family, "binomial"))

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    cv <- do.call("lmerControl", control)

    ## Should difficulties be modeled as random effects?

    ranDiff <- ".item" %in% all.vars(formula)
    ## check for a binary matrix response
    stopifnot(is.matrix(Y), is.numeric(Y),
              all(unique(as.vector(Y)) %in% c(0, 1, NA)))
    nr <- nrow(Y)
    nc <- ncol(Y)
    ## expand model frame etc according to the items
    ind <- rep.int(1:nr, nc)
    mf <- fr$mf[ind, ]
    X <- fr$X[ind, ]
    if (ranDiff) mf$.item <- gl(nc, nr)
    else X <- cbind(X, -contr.sum(nc)[rep(1:nc,each = nr),])
    Y <- as.vector(Y)
    offset <- fr$offset[ind]

    mf$.subj <- gl(nr, 1, nr*nc)
    form <- substitute(Y ~ base + (1|.subj), list(base = formula[[3]]))
    ## establish factor list and Ztl
    FL <- lmerFactorList(form, mf, 0L, 0L)
    fl <- FL$fl
    ## initial fit of a glm to the fixed-effects only.
    glmFit <- glm.fit(X, Y, weights = fl$weights[ind],
                      offset = offset, family = family,
                      intercept = attr(fr$mt, "intercept") > 0)
    Y <- as.double(glmFit$y)
    glmFit

##     ## extract some of the components of glmFit
##     ## weights could have changed
##     weights <- glmFit$prior.weights
##     eta <- glmFit$linear.predictors
##     linkinv <- quote(family$linkinv(eta))
##     mu.eta <- quote(family$mu.eta(eta))
##     mu <- family$linkinv(eta)
##     variance <- quote(family$variance(mu))
##     dev.resids <- quote(family$dev.resids(Y, mu, weights))
##     doLMEopt <- quote(LMEopt(x = mer, value = cv))
##     mer@devComp[8] <- -1
##     mer@status["glmm"] <- as.integer(2) # always use Laplace
##     GSpt <- .Call(glmer_init, environment(), fltype)
##     PQLpars <- c(coef(glmFit), .Call(mer_coef, mer, 2))
##     fixInd <- seq(ncol(X))
##     ## pars[fixInd] == beta, pars[-fixInd] == theta
##     ## indicator of constrained parameters
##     const <- c(rep(FALSE, length(fixInd)),
##                unlist(lapply(mer@nc[seq(along = fl)],
##                              function(k) 1:((k*(k+1))/2) <= k)
##                       ))
## ##     devLaplace <- function(pars) .Call(glmer_devLaplace, pars, GSpt)
##     devLaplace <- NULL
##     rel.tol <- abs(0.01/devLaplace(PQLpars))
##     cat(paste("relative tolerance set to", rel.tol, "\n"))

##     optimRes <- nlminb(PQLpars, devLaplace,
##                        lower = ifelse(const, 5e-10, -Inf),
##                        control = list(trace = cv$msVerbose,
##                        iter.max = cv$msMaxIter,
##                        rel.tol = rel.tol))
##     .Call(glmer_finalize, GSpt)
##     new("glmer",
##         new("lmer", mer,
##             frame =  data.frame(),
##             terms = fr$mt, call = match.call()),
##         weights = weights,
##         family=family)
}
