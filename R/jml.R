asBinaryMatrix <- function(dat)
{
    ## Convert individual columns of a data frame
    asBinaryNumeric <- function(x) {
        if (is.logical(x)) return(as.numeric(x))
        if (is.factor(x)) {
            if (length(levels(x)) != 2)
                stop("factors in the data frame must have exactly 2 levels")
            return(as.numeric(x) - 1L)
        }
    }

    ## Massage the data into a numeric, binary matrix
    if (is.data.frame(dat))
        dat <- do.call("cbind", lapply(dat, asBinaryNumeric))
    dat <- as.matrix(dat)
    storage.mode(dat) <- "double"
    dat
}

sparseRasch <- function(dat, verbose = -1L)
{
    dat <- asBinaryMatrix(dat)
    m <- nrow(dat)
    q <- ncol(dat)
    n <- m * q
    y <- as.vector(dat)
    stopifnot(length(unique(y)) == 2, min(y) == 0, max(y) == 1)

    ## Generate the transpose of the model matrix using the
    ## contr.treatment scheme.  Ability parameters come first then
    ## difficulties then the intercept, which is the logit of the
    ## probability of a correct response by the first subject on the
    ## first question.
    MM <- rBind(as(gl(m, 1, n), "sparseMatrix"),
                as(gl(q, m), "sparseMatrix")[-1, ])
    MM@Dimnames <- vector("list", 2)
    p <- nrow(MM)
    dd <- VecFromNames(dimsNames, "integer")
    dd["n"] <- n
    dd["p"] <- p
    dd["q"] <- q
    dd["lTyp"] <- 1L
    dd["vTyp"] <- 2L
    ans <- new("sparseRasch",
               dims = dd,
               Zt = MM,
               y = as.double(y),
               deviance = VecFromNames(devNames, "numeric"),
               offset = numeric(0),
               L = .Call(mer_create_L, MM),
               fixef = numeric(p),
               mu = numeric(n),
               muEta = numeric(n),
               pWt = numeric(0),
               resid = numeric(n),
               sqrtrWt = numeric(n),
               var = numeric(n))
    .Call(spR_optimize, ans, -1L)
    ans
##     for (i in c(0, seq_len(maxirls))) { # IRLS iterations
##         eta <- crossprod(MM, beta)@x
##         beta_old <- beta
##         mu <- .Call("logit_linkinv", eta, PACKAGE = "stats")
##         swts <- sqrt(1/(mu * (1 - mu)))
##         wtres <- swts * (y - mu)
##         M1@x <- swts * .Call("logit_mu_eta", eta, PACKAGE = "stats")
##         L <- Cholesky(tcrossprod(M1), perm = FALSE, LDL = FALSE)
##         inc1 <- solve(L, wtres, "L")
##         crit <- sqrt((inc1 %*% inc1)/(wtres %*% wtres))
##         beta <- beta + solve(L, inc1, "Lt")
##         if (crit < tol) break
##     }
##     MM@x[] <- 1
##     list(beta = beta, Xt = MM, mu = mu, L = L)
}

asBinaryFrame <- function(dat)
{
#    dat <- asBinaryMatrix(dat)
    m <- nrow(dat)
    q <- ncol(dat)
    y <- as.vector(dat)
    stopifnot(length(unique(y)) == 2, min(y) == 0, max(y) == 1)
    data.frame(y = y, .subj = gl(m, 1, length = m * q),
               .item = gl(q, m))
#    lmer(y ~ (1|.subj) + (1|.item), df, family = binomial, verbose = verbose)
}

