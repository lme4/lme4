## Shared helper functions for the log-link-gaussian mini parameter-recovery
## survey (GH #936: matching up gaussian(link="log") results across
## packages on nlme::Rail, https://github.com/lme4/lme4/issues/936).
##
## Much more limited than ../paramsurvey/: one dataset (Rail), one family
## (gaussian(link="log")), a single scalar random intercept -- so no need
## for the generic multi-term formula/cnms machinery ../paramsurvey/toolkit.R
## has. Formulas are hardcoded per-method below (mgcv needs its own
## s(Rail, bs="re") syntax anyway, so there's no real unification to be had).
##
## Five methods: glmmTMB (reference), joint-phi (R-level devfun, phi as a
## first-class outer bobyqa parameter -- see ../paramsurvey/toolkit.R,
## §7/§11 of ../README_Gamma_GLMMs.md for the general idea), glmer, mgcv::gam,
## MASS::glmmPQL. No PIRLS/moment-phi R-level devfun arm (that's the same
## algorithm as glmer itself, just reimplemented in R -- not useful here).

suppressMessages({
  library(lme4)
  library(glmmTMB)
  library(mgcv)
  library(MASS)
  library(nlme)
  library(minqa)
  library(reformulas)  # nobars()
})
## RTMB backend doesn't implement Gamma yet (as of this glmmTMB dev
## version); force the legacy backend so gaussian and Gamma runs are on
## the same numerical footing (matches 08_prep_crossed.R's toggle)
glmmTMB:::useRTMB(FALSE)

## ---- standardized per-fit result skeleton ----
## sd/sigma (not sd1/sd2/corr/phi as in ../paramsurvey/toolkit.R): this
## example has exactly one scalar random intercept and gaussian's natural
## dispersion-like quantity is sigma itself (residual SD) -- matching how
## GH #936 itself frames the cross-package comparison ("glmer ... reports a
## different sigma").
emptyResult <- function(i, status, msg, time_sec = NA_real_,
                         beta = NA_real_, sd = NA_real_, sigma = NA_real_,
                         negll = NA_real_, singular = NA) {
  list(i = i, status = status, msg = msg, time_sec = time_sec,
       beta = unname(beta), sd = unname(sd), sigma = unname(sigma),
       negll = negll, singular = singular)
}

resultsToDF <- function(results) {
  data.frame(
    i = vapply(results, `[[`, integer(1), "i"),
    status = vapply(results, `[[`, character(1), "status"),
    singular = vapply(results, function(x) isTRUE(x$singular), logical(1)),
    msg = vapply(results, `[[`, character(1), "msg"),
    time_sec = vapply(results, `[[`, numeric(1), "time_sec"),
    beta = vapply(results, `[[`, numeric(1), "beta"),
    sd = vapply(results, `[[`, numeric(1), "sd"),
    sigma = vapply(results, `[[`, numeric(1), "sigma"),
    negll = vapply(results, `[[`, numeric(1), "negll")
  )
}

## ---- common warning-capturing wrapper ----
withWarnings <- function(expr) {
  warn_msgs <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  list(val = val, warn_msgs = warn_msgs)
}

## ---- fit wrappers, one per method ----

fit_glmmTMB_one <- function(i, dat) {
  t0 <- Sys.time()
  r <- withWarnings(glmmTMB(travel ~ 1 + (1 | Rail), data = dat,
                             family = gaussian(link = "log")))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(r$val, "error")) return(emptyResult(i, "error", conditionMessage(r$val), time_sec))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(performance::check_singularity(fit), error = function(e) NA)
  vc <- attr(VarCorr(fit)$cond$Rail, "stddev")
  negll <- tryCatch({
    ll <- as.numeric(logLik(fit))
    if (is.na(ll)) 2 * fit$obj$fn(fit$fit$par) else -2 * ll
  }, error = function(e) NA_real_)
  emptyResult(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
              beta = fixef(fit)$cond[["(Intercept)"]], sd = vc, sigma = sigma(fit),
              negll = negll, singular = singular)
}

## ---- joint-phi devfun, gaussian-specific ----
## Adapted from ../paramsurvey/toolkit.R's make_joint_phi_devfun (Gamma):
## the PIRLS mode-finding loop is already family-generic (dispatches via
## fam$linkinv/variance/mu.eta/dev.resids), so the only Gamma-specific
## part -- the final density-based fit term -- needs swapping from dgamma
## to dnorm (phi = sigma^2 directly, no shape/scale reparameterization
## needed for gaussian). Kept as its own self-contained copy here rather
## than editing the shared ../paramsurvey/toolkit.R, so as not to disturb
## that survey's own (already-run, already-reported) results.
make_joint_phi_devfun_gaussian <- function(form, data, family = gaussian(link = "log")) {
  gm <- glFormula(form, data = data, family = family)
  X <- gm$X
  y <- gm$fr[[1]]
  Zt <- gm$reTrms$Zt
  Lambdat <- gm$reTrms$Lambdat
  Lind <- gm$reTrms$Lind
  thfun <- function(theta) theta[Lind]
  weights <- rep(1, nrow(X))
  offset <- numeric(nrow(X))
  n <- nrow(X); p <- ncol(X); q <- nrow(Zt)
  nth <- length(gm$reTrms$theta)

  fam <- if (is.function(family)) family() else family
  linkinv <- fam$linkinv; variance <- fam$variance; muEta <- fam$mu.eta
  sqDevResid <- fam$dev.resids

  devfun <- function(params) {
    theta <- params[1:nth]
    logphi <- params[nth + 1]
    beta <- params[(nth + 2):length(params)]
    phi <- exp(logphi)

    Lambdat@x[] <- thfun(theta)
    LtZt <- Lambdat %*% Zt
    offb <- offset + as.vector(X %*% beta)
    eta <- numeric(n); mu <- numeric(n)
    updatemu <- function(uu) {
      eta[] <<- offb + as.vector(crossprod(LtZt, uu))
      mu[] <<- linkinv(eta)
      sum(sqDevResid(y, mu, weights)) / phi + sum(uu^2)
    }
    u <- numeric(q)
    olducden <- updatemu(u)
    L <- Matrix::Cholesky(Matrix::tcrossprod(LtZt), perm = FALSE, LDL = FALSE, Imult = 1)
    cvgd <- FALSE
    for (i in 1:60) {
      Whalf <- Matrix::Diagonal(x = sqrt(weights / (phi * variance(mu))))
      LtZtMWhalf <- LtZt %*% (Matrix::Diagonal(x = muEta(eta)) %*% Whalf)
      L <- Matrix::update(L, LtZtMWhalf, 1)
      wtres <- Whalf %*% (y - mu)
      delu <- as.vector(Matrix::solve(L, LtZtMWhalf %*% wtres - u))
      ucden <- updatemu(u + delu)
      if (abs((olducden - ucden) / ucden) < 1e-8) { cvgd <- TRUE; break }
      if (ucden > olducden) {
        for (j in 1:10) {
          ucden <- updatemu(u + (delu <- delu / 2))
          if (ucden <= olducden) break
        }
        if (ucden > olducden) break
      }
      olducden <- ucden
      u <- u + delu
    }
    if (!cvgd) return(1e10)

    ldL2 <- 2 * determinant(L, logarithm = TRUE, sqrt = TRUE)$modulus
    attributes(ldL2) <- NULL

    ## gaussian: phi IS the variance (sigma^2), no shape/scale
    ## reparameterization needed (unlike Gamma's dgamma(shape=1/phi,
    ## scale=mu*phi) -- see ../paramsurvey/toolkit.R)
    fit_term <- -2 * sum(weights * dnorm(y, mean = mu, sd = sqrt(phi), log = TRUE))
    val <- fit_term + sum(u^2) + ldL2
    if (!is.finite(val)) return(1e10)
    val
  }

  glm0 <- glm(nobars(form), data = data, family = family)
  beta_start <- unname(coef(glm0))
  logphi_start <- log(summary(glm0)$dispersion)
  theta_start <- gm$reTrms$theta

  lower_theta <- ifelse(is.finite(gm$reTrms$lower), gm$reTrms$lower, -10)
  upper_theta <- rep(10, nth)

  list(devfun = devfun, nth = nth, p = p, cnms = gm$reTrms$cnms,
       start = c(theta_start, logphi_start, beta_start),
       lower = c(lower_theta, log(1e-4), rep(-20, p)),
       upper = c(upper_theta, log(1e4), rep(20, p)),
       beta_names = colnames(X))
}

fit_jointphi_one <- function(i, dat) {
  t0 <- Sys.time()
  jp <- make_joint_phi_devfun_gaussian(travel ~ 1 + (1 | Rail), dat)
  r <- withWarnings(
    bobyqa(jp$start, jp$devfun, lower = jp$lower, upper = jp$upper,
           control = list(maxfun = 8000))
  )
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(r$val, "error")) return(emptyResult(i, "error", conditionMessage(r$val), time_sec))
  opt <- r$val
  if (!is.finite(opt$fval) || opt$fval >= 1e10) {
    return(emptyResult(i, "degenerate", "hit sentinel / non-finite", time_sec))
  }
  status <- if (opt$ierr != 0) "warning" else if (length(r$warn_msgs) > 0) "warning" else "clean"
  theta <- opt$par[1:jp$nth]
  phi <- exp(opt$par[jp$nth + 1])
  beta <- opt$par[(jp$nth + 2):length(opt$par)]
  sd_re <- unname(theta[1])  # single scalar RE: theta IS the RE sd directly
  singular <- sd_re < 1e-4
  msg <- paste(c(if (opt$ierr != 0) paste("ierr =", opt$ierr), r$warn_msgs), collapse = "; ")
  emptyResult(i, status, msg, time_sec,
              beta = beta[1], sd = sd_re, sigma = sqrt(phi),
              negll = opt$fval, singular = singular)
}

fit_glmer_one <- function(i, dat) {
  t0 <- Sys.time()
  r <- withWarnings(glmer(travel ~ 1 + (1 | Rail), data = dat,
                           family = gaussian(link = "log")))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(r$val, "error")) return(emptyResult(i, "error", conditionMessage(r$val), time_sec))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(isSingular(fit), error = function(e) NA)
  vc <- attr(VarCorr(fit)$Rail, "stddev")
  negll <- tryCatch(-2 * as.numeric(logLik(fit)), error = function(e) NA_real_)
  emptyResult(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
              beta = fixef(fit)[["(Intercept)"]], sd = vc, sigma = sigma(fit),
              negll = negll, singular = singular)
}

fit_mgcv_one <- function(i, dat) {
  t0 <- Sys.time()
  r <- withWarnings(gam(travel ~ 1 + s(Rail, bs = "re"), method = "ML",
                         data = dat, family = gaussian(link = "log")))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(r$val, "error")) return(emptyResult(i, "error", conditionMessage(r$val), time_sec))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  vcomp <- gam.vcomp(fit, conf.lev = 0.95)$vc
  sd_re <- unname(vcomp["s(Rail)", "std.dev"])
  sigma_est <- unname(vcomp["scale", "std.dev"])
  singular <- sd_re < 1e-4
  ## fit$gcv.ubre is mgcv's own ML criterion (negative log marginal
  ## likelihood) for method="ML" fits -- matches glmmTMB's -logLik closely
  ## (verified: 64.586 vs glmmTMB's 64.6 on the real Rail reference fit).
  ## logLik(fit)/logLik.gam()'s own value is NOT the same thing (it
  ## applies an edf-based correction and is not comparable here).
  negll <- 2 * fit$gcv.ubre
  emptyResult(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
              beta = unname(coef(fit)[["(Intercept)"]]), sd = sd_re, sigma = sigma_est,
              negll = negll, singular = singular)
}

## ---- crossed-RE (Rail1 x Rail2) extension ----
## Two random intercepts instead of one, so sd1/sd2 replace the single sd
## field above. q = total number of random-effect levels
## (nlevels(Rail1) + nlevels(Rail2), after droplevels()) -- a property of
## the *data*, identical across methods fit to the same dataset; carried
## through per-fit mainly as a sanity check that no method silently drops
## levels. glmmPQL is skipped here (crossed grouping needs awkward
## formula tricks in nlme); joint-phi is skipped (its devfun currently
## assumes a single scalar RE term).
emptyResultCrossed <- function(i, status, msg, time_sec = NA_real_,
                                beta = NA_real_, sd1 = NA_real_, sd2 = NA_real_,
                                sigma = NA_real_, negll = NA_real_,
                                singular = NA, q = NA_integer_) {
  list(i = i, status = status, msg = msg, time_sec = time_sec,
       beta = unname(beta), sd1 = unname(sd1), sd2 = unname(sd2),
       sigma = unname(sigma), negll = negll, singular = singular, q = q)
}

resultsToDFCrossed <- function(results) {
  data.frame(
    i = vapply(results, `[[`, integer(1), "i"),
    status = vapply(results, `[[`, character(1), "status"),
    singular = vapply(results, function(x) isTRUE(x$singular), logical(1)),
    msg = vapply(results, `[[`, character(1), "msg"),
    time_sec = vapply(results, `[[`, numeric(1), "time_sec"),
    beta = vapply(results, `[[`, numeric(1), "beta"),
    sd1 = vapply(results, `[[`, numeric(1), "sd1"),
    sd2 = vapply(results, `[[`, numeric(1), "sd2"),
    sigma = vapply(results, `[[`, numeric(1), "sigma"),
    negll = vapply(results, `[[`, numeric(1), "negll"),
    q = vapply(results, function(x) as.integer(x$q), integer(1))
  )
}

qFromData <- function(dat) nlevels(droplevels(dat$Rail1)) + nlevels(droplevels(dat$Rail2))

fit_glmmTMB_crossed_one <- function(i, dat, family = gaussian(link = "log")) {
  t0 <- Sys.time()
  r <- withWarnings(glmmTMB(travel ~ 1 + (1 | Rail1) + (1 | Rail2), data = dat,
                             family = family))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  q <- qFromData(dat)
  if (inherits(r$val, "error")) return(emptyResultCrossed(i, "error", conditionMessage(r$val), time_sec, q = q))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(performance::check_singularity(fit), error = function(e) NA)
  vc <- VarCorr(fit)$cond
  sd1 <- attr(vc$Rail1, "stddev")[["(Intercept)"]]
  sd2 <- attr(vc$Rail2, "stddev")[["(Intercept)"]]
  negll <- tryCatch({
    ll <- as.numeric(logLik(fit))
    if (is.na(ll)) 2 * fit$obj$fn(fit$fit$par) else -2 * ll
  }, error = function(e) NA_real_)
  emptyResultCrossed(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
                      beta = fixef(fit)$cond[["(Intercept)"]], sd1 = sd1, sd2 = sd2,
                      sigma = sigma(fit), negll = negll, singular = singular, q = q)
}

fit_glmer_crossed_one <- function(i, dat, family = gaussian(link = "log")) {
  t0 <- Sys.time()
  r <- withWarnings(glmer(travel ~ 1 + (1 | Rail1) + (1 | Rail2), data = dat,
                           family = family))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  q <- qFromData(dat)
  if (inherits(r$val, "error")) return(emptyResultCrossed(i, "error", conditionMessage(r$val), time_sec, q = q))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(isSingular(fit), error = function(e) NA)
  vc <- VarCorr(fit)
  sd1 <- attr(vc$Rail1, "stddev")[["(Intercept)"]]
  sd2 <- attr(vc$Rail2, "stddev")[["(Intercept)"]]
  negll <- tryCatch(-2 * as.numeric(logLik(fit)), error = function(e) NA_real_)
  emptyResultCrossed(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
                      beta = fixef(fit)[["(Intercept)"]], sd1 = sd1, sd2 = sd2,
                      sigma = sigma(fit), negll = negll, singular = singular, q = q)
}

fit_mgcv_crossed_one <- function(i, dat, family = gaussian(link = "log")) {
  t0 <- Sys.time()
  r <- withWarnings(gam(travel ~ 1 + s(Rail1, bs = "re") + s(Rail2, bs = "re"),
                         method = "ML", data = dat, family = family))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  q <- qFromData(dat)
  if (inherits(r$val, "error")) return(emptyResultCrossed(i, "error", conditionMessage(r$val), time_sec, q = q))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  vcomp <- gam.vcomp(fit, conf.lev = 0.95)$vc
  sd1 <- unname(vcomp["s(Rail1)", "std.dev"])
  sd2 <- unname(vcomp["s(Rail2)", "std.dev"])
  sigma_est <- unname(vcomp["scale", "std.dev"])
  singular <- (sd1 < 1e-4) || (sd2 < 1e-4)
  ## see fit_mgcv_one's comment above re: gcv.ubre as the ML criterion
  negll <- 2 * fit$gcv.ubre
  emptyResultCrossed(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
                      beta = unname(coef(fit)[["(Intercept)"]]), sd1 = sd1, sd2 = sd2,
                      sigma = sigma_est, negll = negll, singular = singular, q = q)
}

## ---- "oneway" control case: all levels of Rail1, one level of Rail2 ----
## A genuine one-way design (no crossing at all) subsampled from the same
## crossed simulation -- q is unambiguously nlevels(Rail1), so this is the
## control against which the structured/random crossed-design q_eff
## puzzle can be compared. Uses emptyResultCrossed()'s sd1/sd2/q skeleton
## (sd2 left NA -- there's no second RE term in this model) so its results
## slot directly into the same combined analysis as the crossed cases.
fit_glmmTMB_oneway_one <- function(i, dat, family = gaussian(link = "log")) {
  t0 <- Sys.time()
  r <- withWarnings(glmmTMB(travel ~ 1 + (1 | Rail1), data = dat,
                             family = family))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  q <- nlevels(droplevels(dat$Rail1))
  if (inherits(r$val, "error")) return(emptyResultCrossed(i, "error", conditionMessage(r$val), time_sec, q = q))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(performance::check_singularity(fit), error = function(e) NA)
  sd1 <- attr(VarCorr(fit)$cond$Rail1, "stddev")[["(Intercept)"]]
  negll <- tryCatch({
    ll <- as.numeric(logLik(fit))
    if (is.na(ll)) 2 * fit$obj$fn(fit$fit$par) else -2 * ll
  }, error = function(e) NA_real_)
  emptyResultCrossed(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
                      beta = fixef(fit)$cond[["(Intercept)"]], sd1 = sd1, sd2 = NA_real_,
                      sigma = sigma(fit), negll = negll, singular = singular, q = q)
}

fit_glmer_oneway_one <- function(i, dat, family = gaussian(link = "log")) {
  t0 <- Sys.time()
  r <- withWarnings(glmer(travel ~ 1 + (1 | Rail1), data = dat,
                           family = family))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  q <- nlevels(droplevels(dat$Rail1))
  if (inherits(r$val, "error")) return(emptyResultCrossed(i, "error", conditionMessage(r$val), time_sec, q = q))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(isSingular(fit), error = function(e) NA)
  sd1 <- attr(VarCorr(fit)$Rail1, "stddev")[["(Intercept)"]]
  negll <- tryCatch(-2 * as.numeric(logLik(fit)), error = function(e) NA_real_)
  emptyResultCrossed(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
                      beta = fixef(fit)[["(Intercept)"]], sd1 = sd1, sd2 = NA_real_,
                      sigma = sigma(fit), negll = negll, singular = singular, q = q)
}

fit_mgcv_oneway_one <- function(i, dat, family = gaussian(link = "log")) {
  t0 <- Sys.time()
  r <- withWarnings(gam(travel ~ 1 + s(Rail1, bs = "re"), method = "ML",
                         data = dat, family = family))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  q <- nlevels(droplevels(dat$Rail1))
  if (inherits(r$val, "error")) return(emptyResultCrossed(i, "error", conditionMessage(r$val), time_sec, q = q))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  vcomp <- gam.vcomp(fit, conf.lev = 0.95)$vc
  sd1 <- unname(vcomp["s(Rail1)", "std.dev"])
  sigma_est <- unname(vcomp["scale", "std.dev"])
  singular <- sd1 < 1e-4
  negll <- 2 * fit$gcv.ubre
  emptyResultCrossed(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
                      beta = unname(coef(fit)[["(Intercept)"]]), sd1 = sd1, sd2 = NA_real_,
                      sigma = sigma_est, negll = negll, singular = singular, q = q)
}

fit_pql_one <- function(i, dat) {
  t0 <- Sys.time()
  r <- withWarnings(glmmPQL(travel ~ 1, random = ~1 | Rail,
                             family = gaussian(link = "log"), data = dat,
                             verbose = FALSE))
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(r$val, "error")) return(emptyResult(i, "error", conditionMessage(r$val), time_sec))
  fit <- r$val
  status <- if (length(r$warn_msgs) > 0) "warning" else "clean"
  sd_re <- as.numeric(nlme::VarCorr(fit)["(Intercept)", "StdDev"])
  singular <- sd_re < 1e-4
  ## PQL is a penalized quasi-likelihood method, not full ML -- it has no
  ## meaningful marginal logLik comparable to the other three methods'
  ## (glmmPQL's own logLik()/AIC()/BIC() are NA, consistent with this).
  emptyResult(i, status, paste(r$warn_msgs, collapse = "; "), time_sec,
              beta = nlme::fixef(fit)[["(Intercept)"]], sd = sd_re, sigma = fit$sigma,
              negll = NA_real_, singular = singular)
}
