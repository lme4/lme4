## Simulation study: which moment-based dispersion denominator for
## free-dispersion GLMMs (glmerControl(disp_dof_correction = "n", "n-p",
## "n-q")) gives the least biased / lowest-MSE estimates, compared with
## glmmTMB ML and REML? See misc/TODO_AGQ.md ("the deviance uses two
## different phi values").

suppressPackageStartupMessages({
    library(lme4)
    library(glmmTMB)
})

## Gamma(log) GLMM, y ~ x + (1|g). phi = 1/shape is the true dispersion.
simdata <- function(ngrp = 20, nper = 20, beta = c(1, 0.5), sd_re = 0.5,
                    shape = 0.5) {
    n <- ngrp * nper
    dd <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ngrp), each = nper)))
    b <- rnorm(ngrp, 0, sd_re)
    mu <- exp(beta[1] + beta[2] * dd$x + b[dd$g])
    dd$y <- rgamma(n, shape = shape, scale = mu / shape)
    dd
}

## one row per method: estimates, logLik, timing, warnings
fit_all <- function(dd, form = y ~ x + (1 | g)) {
    res <- list()
    wrap <- function(method, expr) {
        w <- character(0)
        t0 <- proc.time()[["elapsed"]]
        fit <- withCallingHandlers(
            tryCatch(expr, error = function(e) e),
            warning = function(cnd) {
                w <<- c(w, conditionMessage(cnd))
                invokeRestart("muffleWarning")
            })
        time <- proc.time()[["elapsed"]] - t0
        list(method = method, fit = fit, time = time,
             warn = paste(unique(w), collapse = " | "))
    }
    glmer_dof <- function(type, start = NULL)
        glmer(form, dd, Gamma(link = "log"), start = start,
              control = glmerControl(disp_dof_correction = type))
    ## warm start the corrected fits from the uncorrected one
    r_n <- wrap("glmer_n", glmer_dof("n"))
    start <- if (inherits(r_n$fit, "merMod"))
                 list(theta = getME(r_n$fit, "theta"), fixef = fixef(r_n$fit))
    res <- list(r_n,
                wrap("glmer_n-p", glmer_dof("n-p", start)),
                wrap("glmer_n-q", glmer_dof("n-q", start)),
                wrap("glmmTMB_ML", glmmTMB(form, dd, Gamma(link = "log"),
                                           control = glmmTMBControl(parallel = 1))),
                wrap("glmmTMB_REML", glmmTMB(form, dd, Gamma(link = "log"), REML = TRUE,
                                             control = glmmTMBControl(parallel = 1))))
    do.call(rbind, lapply(res, extract_row))
}

extract_row <- function(r) {
    out <- data.frame(method = r$method, beta0 = NA_real_, beta1 = NA_real_,
                      sd_re = NA_real_, phi = NA_real_, logLik = NA_real_,
                      singular = NA, time = r$time, warn = r$warn, error = "")
    f <- r$fit
    if (inherits(f, "error")) {
        out$error <- conditionMessage(f)
    } else if (inherits(f, "merMod")) {
        out[c("beta0", "beta1")] <- fixef(f)
        out$sd_re <- getME(f, "theta")[[1]]     # Gamma: theta is the RE sd
        out$phi <- sigma(f)^2
        out$logLik <- c(logLik(f))
        out$singular <- isSingular(f)
        out <- cbind(out, phi_alt(f))
    } else {
        out[c("beta0", "beta1")] <- fixef(f)$cond
        out$sd_re <- sqrt(VarCorr(f)$cond$g[1, 1])
        out$phi <- sigma(f)^2                   # glmmTMB Gamma: sigma^2 = 1/shape
        out$logLik <- c(logLik(f))
        out$singular <- out$sd_re < 1e-4
    }
    if (is.null(out$phi_digamma))
        out[c("phi_digamma", "phi_digamma_q", "phi_pearson_p", "phi_pearson_q")] <- NA_real_
    out
}

## Reporting-only alternatives to sigma(fit)^2 for a Gamma glmer fit, computed
## post hoc at the fit's conditional means (so theta/beta are unchanged):
##  - phi_digamma: conditional ML given mu, 1/nu where
##    log(nu) - digamma(nu) = D/(2n) (as in MASS::gamma.shape)
##  - phi_digamma_q: the same with D/(2(n - q)), an ad hoc df correction
##  - phi_pearson_p, phi_pearson_q: Pearson X^2/(n - p) and X^2/(n - q),
##    q = rank of [X, Z]
phi_alt <- function(f) {
    y <- getME(f, "y"); mu <- fitted(f); n <- length(y)
    D <- sum(Gamma()$dev.resids(y, mu, 1))
    nu_root <- function(m) uniroot(function(nu) log(nu) - digamma(nu) - D / (2 * m),
                                   c(1e-8, 1e8), tol = 1e-12)$root
    X <- getME(f, "X")
    q <- rankMatrix(cbind(X, getME(f, "Z")), method = "qr")[1]
    X2 <- sum((y - mu)^2 / mu^2)
    data.frame(phi_digamma = 1 / nu_root(n), phi_digamma_q = 1 / nu_root(n - q),
               phi_pearson_p = X2 / (n - ncol(X)),
               phi_pearson_q = X2 / (n - q))
}

## one replicate, reproducible regardless of how replicates are scheduled
run_one <- function(i, seed0 = 1000, ...) {
    set.seed(seed0 + i)
    dd <- simdata(...)
    cbind(sim = i, fit_all(dd))
}
