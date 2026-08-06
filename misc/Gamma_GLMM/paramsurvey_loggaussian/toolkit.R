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
## Four methods (not six -- no joint-phi/PIRLS-phi R-level devfun arms this
## time, per this session's request): glmmTMB (reference), glmer, mgcv::gam,
## MASS::glmmPQL.

suppressMessages({
  library(lme4)
  library(glmmTMB)
  library(mgcv)
  library(MASS)
  library(nlme)
})

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
