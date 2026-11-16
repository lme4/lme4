## Adaptive Gauss-Hermite quadrature (nAGQ > 1) for a single vector-valued
## random-effects term, e.g. (1 + x | g): glmerAGQ() in src/external.cpp,
## AGQrule() in R/GHrule.R.

## independent reference for -2 log L: per-level adaptive product
## Gauss-Hermite quadrature on the b scale, with modes from optim() and
## the observed Hessian from optimHess()
refDev <- function(theta, beta, X, Z, y, g, ll_fun, nq = 20, offset = 0) {
    k <- ncol(Z)
    L <- matrix(0, k, k); L[lower.tri(L, TRUE)] <- theta
    Sig <- tcrossprod(L)
    Sinv <- solve(Sig); ldS <- c(determinant(Sig)$modulus)
    gh <- GHrule(nq)
    zz <- as.matrix(expand.grid(rep(list(gh[, "z"]), k)))
    ww <- apply(as.matrix(expand.grid(rep(list(gh[, "w"]), k))), 1, prod)
    eta0 <- drop(X %*% beta) + offset
    tot <- 0
    for (lev in levels(g)) {
        i <- g == lev
        f <- function(b) -(sum(ll_fun(y[i], eta0[i] + Z[i, , drop = FALSE] %*% b)) -
                           0.5 * sum(b * (Sinv %*% b)) - 0.5 * ldS - k/2 * log(2*pi))
        o <- optim(rep(0, k), f, method = "BFGS", control = list(reltol = 1e-12))
        R <- chol(optimHess(o$par, f))
        bb <- t(o$par + backsolve(R, diag(k)) %*% t(zz))
        lv <- -apply(bb, 1, f) - rowSums(dnorm(zz, log = TRUE))
        m <- max(lv)
        tot <- tot + m + log(sum(ww * exp(lv - m))) - sum(log(diag(R)))
    }
    -2 * tot
}

set.seed(101)
ng <- 20; nper <- 8
dd <- data.frame(g = gl(ng, nper), x = rnorm(ng * nper))
bb <- cbind(rnorm(ng), rnorm(ng, sd = 0.7))
eta <- -0.5 + 0.7 * dd$x + bb[dd$g, 1] + bb[dd$g, 2] * dd$x
dd$yb <- rbinom(nrow(dd), 1, plogis(eta))
dd$yp <- rpois(nrow(dd), exp(eta))
theta0 <- c(0.9, 0.2, 0.6); beta0 <- c(-0.4, 0.6)
XX <- cbind(1, dd$x)

test_that("AGQrule() reduces to GHrule() and integrates low-order moments", {
    expect_identical(lme4:::AGQrule(7), GHrule(7))
    sp1 <- lme4:::AGQrule(5, 1, "sparse")
    expect_equal(unname(sp1[order(sp1[, 1]), ]), unname(GHrule(5)))
    for (type in c("product", "sparse")) {
        r <- lme4:::AGQrule(5, 3, type)
        expect_equal(ncol(r), 5L)
        w <- r[, "w"]; z <- r[, 1:3]
        expect_equal(sum(w), 1)
        expect_equal(unname(colSums(w * z^2)), rep(1, 3))
        expect_equal(unname(colSums(w * z^4)), rep(3, 3))
        expect_equal(sum(w * z[, 1]^2 * z[, 2]^2), 1)
        expect_equal(r[, "ldnorm"], rowSums(dnorm(z, log = TRUE)))
    }
    expect_equal(nrow(lme4:::AGQrule(5, 2)), 25L)
    expect_lt(nrow(lme4:::AGQrule(5, 3, prune = 1e-4)), 125L)
})

test_that("vector-valued AGQ deviance converges to an independent reference", {
    cases <- list(
        binomial = list(f = yb ~ x + (x|g), fam = binomial, y = dd$yb,
                        ll = function(y, eta) dbinom(y, 1, plogis(eta), log = TRUE)),
        poisson  = list(f = yp ~ x + (x|g), fam = poisson, y = dd$yp,
                        ll = function(y, eta) dpois(y, exp(eta), log = TRUE)))
    for (cc in cases) {
        ref <- refDev(theta0, beta0, XX, XX, cc$y, dd$g, cc$ll)
        dev <- function(n, ...) {
            df <- glmer(cc$f, dd, cc$fam, nAGQ = n, devFunOnly = TRUE, ...)
            df(c(theta0, beta0))
        }
        err <- abs(sapply(c(3, 15), dev) - ref)
        expect_lt(err[2], 1e-4)
        expect_lt(err[2], err[1])
        expect_lt(abs(dev(9, control = glmerControl(AGQ_grid = "sparse")) - ref), 0.01)
    }
})

## do.call() so that weights/offset vectors passed through ... reach
## model.frame() as values (glmer() evaluates them non-standardly)
devAt <- function(formula, data, family, nAGQ, ...)
    do.call(glmer, list(formula, data, family, nAGQ = nAGQ, devFunOnly = TRUE,
                        ...))(c(theta0, beta0))

test_that("vector-valued AGQ handles binomial sizes and prior weights", {
    set.seed(102)
    db <- dd
    db$n <- sample(1:6, nrow(db), replace = TRUE)
    db$s <- rbinom(nrow(db), db$n, plogis(eta))
    db$p <- db$s / db$n
    long <- db[rep(seq_len(nrow(db)), db$n), ]
    long$y <- unlist(lapply(seq_len(nrow(db)), function(i)
        rep(1:0, c(db$s[i], db$n[i] - db$s[i]))))
    ref <- refDev(theta0, beta0, XX, XX, seq_len(nrow(db)), db$g,
                  function(i, eta) dbinom(db$s[i], db$n[i], plogis(eta), log = TRUE))
    for (nq in c(3, 15)) {
        d_cbind <- devAt(cbind(s, n - s) ~ x + (x|g), db, binomial, nq)
        expect_equal(devAt(p ~ x + (x|g), db, binomial, nq, weights = db$n), d_cbind)
        ## the 0/1 data differ only by the binomial coefficients
        expect_equal(devAt(y ~ x + (x|g), long, binomial, nq) - d_cbind,
                     2 * sum(lchoose(db$n, db$s)))
    }
    expect_lt(abs(d_cbind - ref), 1e-4)
})

test_that("vector-valued AGQ handles offsets", {
    set.seed(103)
    do <- dd
    do$expo <- runif(nrow(do), 0.5, 3)
    do$y <- rpois(nrow(do), do$expo * exp(eta))
    ref <- refDev(theta0, beta0, XX, XX, do$y, do$g,
                  function(y, eta) dpois(y, exp(eta), log = TRUE),
                  offset = log(do$expo))
    d_form <- devAt(y ~ x + offset(log(expo)) + (x|g), do, poisson, 15)
    expect_equal(devAt(y ~ x + (x|g), do, poisson, 15, offset = log(do$expo)), d_form)
    expect_lt(abs(d_form - ref), 1e-4)
    m <- glmer(y ~ x + offset(log(expo)) + (x|g), do, poisson, nAGQ = 5)
    expect_equal(deviance(refit(m)), deviance(m), tolerance = 1e-5)
})

test_that("vector-valued AGQ works with nAGQ0initStep = FALSE", {
    fit <- function(init) glmer(yp ~ x + (x|g), dd, poisson, nAGQ = 5,
                                control = glmerControl(nAGQ0initStep = init))
    m_t <- fit(TRUE); m_f <- fit(FALSE)
    expect_equal(c(getME(m_f, "theta"), fixef(m_f)),
                 c(getME(m_t, "theta"), fixef(m_t)), tolerance = 1e-3)
    expect_equal(deviance(m_f), deviance(m_t), tolerance = 1e-5)
})

test_that("sparse-grid fits warn when the quadrature looks unreliable", {
    df <- glmer(yb ~ x + (x|g), dd, binomial, nAGQ = 5, devFunOnly = TRUE,
                control = glmerControl(AGQ_grid = "sparse"))
    df(c(theta0, beta0))
    corr <- environment(df)$AGQcorr
    expect_length(corr, nlevels(dd$g))
    expect_silent(lme4:::glmerAGQcheck(environment(df)))
    chk <- function(grid, corr)
        lme4:::glmerAGQcheck(list2env(list(AGQ_grid = grid, AGQcorr = corr)))
    expect_warning(chk("sparse", rep(0.5, 10)), "below the Laplace")
    expect_warning(chk("sparse", c(NA, rep(0, 9))), "non-positive")
    expect_silent(chk("product", rep(0.5, 10)))
})

test_that("vector-valued AGQ fits, prints and refits", {
    m1 <- glmer(yb ~ x + (x|g), dd, binomial)
    m5 <- glmer(yb ~ x + (x|g), dd, binomial, nAGQ = 5)
    expect_equal(fixef(m5), fixef(m1), tolerance = 0.1)
    expect_match(summary(m5)$methTitle, "nAGQ = 5, 25 nodes", fixed = TRUE)
    expect_equal(deviance(refit(m5)), deviance(m5), tolerance = 1e-5)
    ms <- glmer(yb ~ x + (x|g), dd, binomial, nAGQ = 5,
                control = glmerControl(AGQ_grid = "sparse"))
    expect_match(summary(ms)$methTitle, "sparse grid", fixed = TRUE)
    expect_equal(deviance(ms), deviance(m5), tolerance = 1e-3)
    expect_equal(deviance(refit(ms)), deviance(ms), tolerance = 1e-5)
    ## scalar fits keep the old title
    m0 <- glmer(yb ~ x + (1|g), dd, binomial, nAGQ = 5)
    expect_match(summary(m0)$methTitle, "nAGQ = 5)", fixed = TRUE)
})

test_that("AGQ checks the random-effects structure and grid size", {
    expect_error(glmer(yb ~ x + (x||g), dd, binomial, nAGQ = 3),
                 "single random-effects term")
    dd$h <- gl(4, 40)
    expect_error(glmer(yb ~ x + (1|g) + (1|h), dd, binomial, nAGQ = 3),
                 "single random-effects term")
    dd$z <- rnorm(nrow(dd))
    expect_warning(glmer(yp ~ x + z + (x + z|g), dd, poisson, nAGQ = 11,
                         devFunOnly = TRUE),
                   "1331 quadrature nodes")
    expect_error(glmerControl(AGQ_prune = 2), "AGQ_prune")
})
