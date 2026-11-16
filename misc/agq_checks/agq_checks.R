## Checks of vector-valued AGQ (glmer(nAGQ > 1) with (x|g)) for cases not
## covered by tests/testthat/test-agq-vector.R: binomial response
## specifications and prior weights, offsets, nAGQ0initStep = FALSE, and
## sparse grids on strongly non-Gaussian integrands.
## See misc/TODO_AGQ.md.

library(lme4)
set.seed(20261115)

## independent reference for -2 log L: per-level adaptive product
## Gauss-Hermite quadrature on the b scale (as in test-agq-vector.R), with
## a per-observation offset. ll_fun(i, eta) returns the log-likelihood
## contributions of observations i at linear predictor eta.
refDev <- function(theta, beta, X, Z, g, ll_fun, offset = 0, nq = 20) {
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
        i <- which(g == lev)
        f <- function(b) -(sum(ll_fun(i, eta0[i] + Z[i, , drop = FALSE] %*% b)) -
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

## do.call() so that weights/offset vectors passed through ... reach
## model.frame() as values (glmer() evaluates them non-standardly)
devAt <- function(formula, data, family, pars, nAGQ, ...)
    do.call(glmer, list(formula, data, family, nAGQ = nAGQ, devFunOnly = TRUE, ...))(pars)

ng <- 20; nper <- 8
dd <- data.frame(g = gl(ng, nper), x = rnorm(ng * nper))
bb <- cbind(rnorm(ng), rnorm(ng, sd = 0.7))
eta <- -0.5 + 0.7 * dd$x + bb[dd$g, 1] + bb[dd$g, 2] * dd$x
XX <- cbind(1, dd$x)
theta0 <- c(0.9, 0.2, 0.6); beta0 <- c(-0.4, 0.6); pars0 <- c(theta0, beta0)
nq_vec <- c(1, 5, 11)

## ---------------------------------------------------------------------
## 1. binomial: cbind(s, n - s) vs s/n + weights vs disaggregated 0/1
## ---------------------------------------------------------------------
dd$n <- sample(1:6, nrow(dd), replace = TRUE)
dd$s <- rbinom(nrow(dd), dd$n, plogis(eta))
dd$p <- dd$s / dd$n
long <- dd[rep(seq_len(nrow(dd)), dd$n), ]
long$yb <- unlist(lapply(seq_len(nrow(dd)), function(i)
    rep(1:0, c(dd$s[i], dd$n[i] - dd$s[i]))))
lchoose_const <- 2 * sum(lchoose(dd$n, dd$s)) # 0/1 data have no lchoose term

ref_binom <- refDev(theta0, beta0, XX, XX, dd$g,
                    function(i, eta) dbinom(dd$s[i], dd$n[i], plogis(eta), log = TRUE))
binom_tab <- t(sapply(nq_vec, function(nq) {
    d_cbind <- devAt(cbind(s, n - s) ~ x + (x|g), dd, binomial, pars0, nq)
    d_wts   <- devAt(p ~ x + (x|g), dd, binomial, pars0, nq, weights = dd$n)
    d_long  <- devAt(yb ~ x + (x|g), long, binomial, pars0, nq)
    c(nAGQ = nq, cbind_minus_wts = d_cbind - d_wts,
      long_minus_cbind_minus_const = d_long - d_cbind - lchoose_const,
      cbind_minus_ref = d_cbind - ref_binom)
}))
cat("\n## binomial specifications (ref =", ref_binom, ")\n")
print(signif(binom_tab, 4))

## fitted models agree too
fits_binom <- list(
    cbind = glmer(cbind(s, n - s) ~ x + (x|g), dd, binomial, nAGQ = 7),
    wts   = glmer(p ~ x + (x|g), dd, binomial, nAGQ = 7, weights = n),
    long  = glmer(yb ~ x + (x|g), long, binomial, nAGQ = 7))
print(sapply(fits_binom, function(m) c(getME(m, "theta"), fixef(m),
                                       dev = deviance(m), logLik = c(logLik(m)))))
## refit() keeps the weights
print(sapply(fits_binom[1:2], function(m) deviance(refit(m)) - deviance(m)))

## Poisson with non-integer prior weights (weighted log-likelihood)
dd$yp <- rpois(nrow(dd), exp(eta))
dd$w <- runif(nrow(dd), 0.5, 2)
ref_poisw <- refDev(theta0, beta0, XX, XX, dd$g,
                    function(i, eta) dd$w[i] * dpois(dd$yp[i], exp(eta), log = TRUE))
## R's poisson()$aic for non-integer y*w uses lgamma, so the brute force
## uses the same weighted log-density; for comparison compute the constant
## offset between the two at each nAGQ, which should be independent of nAGQ
## and match the difference at nAGQ = 1 (Laplace, already established)
poisw_tab <- t(sapply(nq_vec, function(nq)
    c(nAGQ = nq,
      dev_minus_ref = suppressWarnings(devAt(yp ~ x + (x|g), dd, poisson, pars0, nq,
                                              weights = dd$w)) - ref_poisw)))
cat("\n## Poisson, prior weights\n")
print(signif(poisw_tab, 4))

## ---------------------------------------------------------------------
## 2. offsets
## ---------------------------------------------------------------------
dd$expo <- runif(nrow(dd), 0.5, 3)
dd$ypo <- rpois(nrow(dd), dd$expo * exp(eta))
ref_off <- refDev(theta0, beta0, XX, XX, dd$g,
                  function(i, eta) dpois(dd$ypo[i], exp(eta), log = TRUE),
                  offset = log(dd$expo))
dd$off_b <- rnorm(nrow(dd), sd = 0.5)
dd$ybo <- rbinom(nrow(dd), 1, plogis(eta + dd$off_b))
ref_offb <- refDev(theta0, beta0, XX, XX, dd$g,
                   function(i, eta) dbinom(dd$ybo[i], 1, plogis(eta), log = TRUE),
                   offset = dd$off_b)
off_tab <- t(sapply(nq_vec, function(nq) {
    d_form <- devAt(ypo ~ x + offset(log(expo)) + (x|g), dd, poisson, pars0, nq)
    d_arg  <- devAt(ypo ~ x + (x|g), dd, poisson, pars0, nq, offset = log(dd$expo))
    d_bin  <- devAt(ybo ~ x + offset(off_b) + (x|g), dd, binomial, pars0, nq)
    c(nAGQ = nq, form_minus_arg = d_form - d_arg, pois_minus_ref = d_form - ref_off,
      binom_minus_ref = d_bin - ref_offb)
}))
cat("\n## offsets\n")
print(signif(off_tab, 4))
m_off <- glmer(ypo ~ x + offset(log(expo)) + (x|g), dd, poisson, nAGQ = 7)
cat("refit() with offset, deviance change:", deviance(refit(m_off)) - deviance(m_off), "\n")

## ---------------------------------------------------------------------
## 3. nAGQ0initStep = FALSE
## ---------------------------------------------------------------------
dd$yb <- rbinom(nrow(dd), 1, plogis(eta))
init_tab <- do.call(rbind, lapply(c("product", "sparse"), function(grid) {
    do.call(rbind, lapply(list(binomial = list(yb ~ x + (x|g), binomial),
                               poisson  = list(yp ~ x + (x|g), poisson)), function(cc) {
        fit <- function(init) glmer(cc[[1]], dd, cc[[2]], nAGQ = 5,
                                    control = glmerControl(nAGQ0initStep = init,
                                                           AGQ_grid = grid))
        m_t <- fit(TRUE); m_f <- fit(FALSE)
        data.frame(family = cc[[2]]()$family, grid = grid,
                   dev_diff = deviance(m_f) - deviance(m_t),
                   max_par_diff = max(abs(c(getME(m_f, "theta") - getME(m_t, "theta"),
                                            fixef(m_f) - fixef(m_t)))),
                   refit_diff = deviance(refit(m_f)) - deviance(m_f))
    }))
}))
cat("\n## nAGQ0initStep = FALSE vs TRUE (nAGQ = 5)\n")
print(init_tab, digits = 3)

## ---------------------------------------------------------------------
## 4. sparse grids on strongly non-Gaussian integrands
## ---------------------------------------------------------------------
## Binary data, small clusters, large RE sd, k = 2 or 3. For a grid of
## theta values (multiples of the true theta) at the true beta, compare
## the sparse-grid devfun with a high-order product rule, and count
## sentinel values (1e10: some level's quadrature sum <= 0). Then fit
## with the sparse grid, counting sentinel hits during optimization.
sentinel <- 1e10
## small clusters with k = 3 have fewer observations than random effects
ctrl0 <- glmerControl(check.nobs.vs.nRE = "ignore")
simBin <- function(ng, nper, k, sd_re, seed) {
    set.seed(seed)
    d <- data.frame(g = gl(ng, nper), x = rnorm(ng * nper), z = rnorm(ng * nper))
    Z <- cbind(1, d$x, d$z)[, seq_len(k), drop = FALSE]
    b <- matrix(rnorm(ng * k, sd = sd_re), ng, k)
    d$y <- rbinom(nrow(d), 1, plogis(-0.3 + 0.5 * d$x + rowSums(Z * b[d$g, , drop = FALSE])))
    d
}
form_k <- list(`2` = y ~ x + (x|g), `3` = y ~ x + (x + z|g))
nq_ref <- c(`2` = 21, `3` = 11)

sparse_grid <- expand.grid(nper = c(2, 3, 4), k = 2:3, sd_re = c(1, 2.5), rep = 1:2)
sparse_res <- do.call(rbind, lapply(seq_len(nrow(sparse_grid)), function(r) {
    cc <- sparse_grid[r, ]
    d <- simBin(40, cc$nper, cc$k, cc$sd_re, seed = 1000 + r)
    f <- form_k[[as.character(cc$k)]]
    th_true <- diag(cc$k)[lower.tri(diag(cc$k), TRUE)] * cc$sd_re
    mk <- function(nq, grid) suppressWarnings(
        glmer(f, d, binomial, nAGQ = nq, devFunOnly = TRUE,
              control = glmerControl(AGQ_grid = grid, check.nobs.vs.nRE = "ignore")))
    dref <- mk(nq_ref[[as.character(cc$k)]], "product")
    do.call(rbind, lapply(c(3, 5, 7, 9), function(nq) {
        dsp <- mk(nq, "sparse"); dpr <- mk(nq, "product")
        vals <- t(sapply(c(0.25, 0.5, 1, 2, 4), function(s) {
            p <- c(s * th_true, -0.3, 0.5)
            c(ref = dref(p), sparse = dsp(p), product = dpr(p))
        }))
        ok <- vals[, "sparse"] < sentinel
        data.frame(cc, nAGQ = nq,
                   n_sentinel = sum(!ok),
                   max_err_sparse = if (any(ok)) max(abs(vals[ok, "sparse"] - vals[ok, "ref"])) else NA,
                   max_err_product = max(abs(vals[, "product"] - vals[, "ref"])),
                   row.names = NULL)
    }))
}))
cat("\n## sparse vs product rule on the devfun, 5 theta values per row\n")
print(sparse_res, digits = 3)
print(aggregate(cbind(n_sentinel, max_err_sparse, max_err_product) ~ k + nAGQ,
                sparse_res, function(x) mean(x, na.rm = TRUE)), digits = 3)

## full fits with sparse grids: count sentinel hits during optimization
fit_res <- do.call(rbind, lapply(seq_len(nrow(sparse_grid)), function(r) {
    cc <- sparse_grid[r, ]
    d <- simBin(40, cc$nper, cc$k, cc$sd_re, seed = 1000 + r)
    f <- form_k[[as.character(cc$k)]]
    do.call(rbind, lapply(c(3, 5), function(nq) {
        df <- suppressWarnings(glmer(f, d, binomial, nAGQ = nq, devFunOnly = TRUE,
                                     control = glmerControl(AGQ_grid = "sparse", check.nobs.vs.nRE = "ignore")))
        hits <- 0; n_eval <- 0
        df2 <- function(p) { v <- df(p); n_eval <<- n_eval + 1
            if (v >= sentinel) hits <<- hits + 1; v }
        m1 <- glmer(f, d, binomial, control = ctrl0)
        rho <- environment(df)
        opt <- Nelder_Mead(df2, c(getME(m1, "theta"), fixef(m1)),
                           lower = rho$lower, upper = rho$upper)
        m_pr <- suppressWarnings(glmer(f, d, binomial, nAGQ = nq, control = ctrl0))
        data.frame(cc, nAGQ = nq, n_eval = n_eval, sentinel_hits = hits,
                   dev_sparse_fit = opt$fval, dev_product_fit = -2 * c(logLik(m_pr)),
                   row.names = NULL)
    }))
}))
cat("\n## optimization with sparse grids (Nelder-Mead from the Laplace fit)\n")
print(fit_res, digits = 4)

## ---------------------------------------------------------------------
## 5. glmer() fits: accuracy at the estimates, and calibration of the
##    sparse-grid warning in glmerAGQcheck() (R/lmer.R)
## ---------------------------------------------------------------------
## For each fit: -2 log L at the estimates from the fit and from a
## high-order product rule, the per-level AGQ corrections (attribute
## "AGQcorr" of the devfun value) and whether glmer() warned.
calib <- do.call(rbind, lapply(seq_len(nrow(sparse_grid)), function(r) {
    cc <- sparse_grid[r, ]
    d <- simBin(40, cc$nper, cc$k, cc$sd_re, seed = 1000 + r)
    f <- form_k[[as.character(cc$k)]]
    dref <- suppressWarnings(glmer(f, d, binomial, nAGQ = nq_ref[[as.character(cc$k)]],
                                   devFunOnly = TRUE, control = ctrl0))
    do.call(rbind, lapply(c(3, 5, 9), function(nq) do.call(rbind, lapply(c("product", "sparse"), function(grid) {
        ctrl <- glmerControl(AGQ_grid = grid, check.nobs.vs.nRE = "ignore")
        w <- character(0)
        m <- withCallingHandlers(glmer(f, d, binomial, nAGQ = nq, control = ctrl),
                                 warning = function(x) {
                                     w <<- c(w, conditionMessage(x))
                                     invokeRestart("muffleWarning")
                                 })
        p <- c(getME(m, "theta"), fixef(m))
        df <- suppressWarnings(glmer(f, d, binomial, nAGQ = nq, devFunOnly = TRUE, control = ctrl))
        df(p); corr <- environment(df)$AGQcorr
        data.frame(cc, nAGQ = nq, grid = grid,
                   err = -2 * c(logLik(m)) - dref(p),
                   max_corr = max(abs(corr)), mean_corr = mean(corr),
                   warned = any(grepl("sparse-grid quadrature", w)), row.names = NULL)
    }))))
}))
cat("\n## glmer() fits: error in -2 log L at the estimates, AGQ corrections\n")
print(aggregate(cbind(abs(err), max_corr, abs(mean_corr)) ~ grid + nAGQ, calib, max), digits = 3)
cat("fits with |err| > 5 or a warning:\n")
print(calib[abs(calib$err) > 5 | calib$warned, ], digits = 3)

## ---------------------------------------------------------------------
## 5. glmer() fits: accuracy at the estimates, and calibration of the
##    sparse-grid warning in glmerAGQcheck() (R/lmer.R)
## ---------------------------------------------------------------------
## For each fit: -2 log L at the estimates from the fit and from a
## high-order product rule, the per-level AGQ corrections (attribute
## "AGQcorr" of the devfun value) and whether glmer() warned.
calib <- do.call(rbind, lapply(seq_len(nrow(sparse_grid)), function(r) {
    cc <- sparse_grid[r, ]
    d <- simBin(40, cc$nper, cc$k, cc$sd_re, seed = 1000 + r)
    f <- form_k[[as.character(cc$k)]]
    dref <- suppressWarnings(glmer(f, d, binomial, nAGQ = nq_ref[[as.character(cc$k)]],
                                   devFunOnly = TRUE, control = ctrl0))
    do.call(rbind, lapply(c(3, 5, 9), function(nq) do.call(rbind, lapply(c("product", "sparse"), function(grid) {
        ctrl <- glmerControl(AGQ_grid = grid, check.nobs.vs.nRE = "ignore")
        w <- character(0)
        m <- withCallingHandlers(glmer(f, d, binomial, nAGQ = nq, control = ctrl),
                                 warning = function(x) {
                                     w <<- c(w, conditionMessage(x))
                                     invokeRestart("muffleWarning")
                                 })
        p <- c(getME(m, "theta"), fixef(m))
        df <- suppressWarnings(glmer(f, d, binomial, nAGQ = nq, devFunOnly = TRUE, control = ctrl))
        df(p); corr <- environment(df)$AGQcorr
        data.frame(cc, nAGQ = nq, grid = grid,
                   err = -2 * c(logLik(m)) - dref(p),
                   max_corr = max(abs(corr)), mean_corr = mean(corr),
                   warned = any(grepl("sparse-grid quadrature", w)), row.names = NULL)
    }))))
}))
calib$abs_err <- abs(calib$err)
calib$abs_mean_corr <- abs(calib$mean_corr)
cat("\n## glmer() fits: error in -2 log L at the estimates, AGQ corrections (max over fits)\n")
print(aggregate(cbind(abs_err, max_corr, abs_mean_corr) ~ grid + nAGQ, calib, max), digits = 3)
cat("fits with |err| > 5 or a warning:\n")
print(calib[calib$abs_err > 5 | calib$warned, ], digits = 3)
