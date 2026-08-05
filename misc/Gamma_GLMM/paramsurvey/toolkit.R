## Shared helper functions for the parameter-recovery survey (broader
## successor to the multistart diagnostic): for a handful of real datasets,
## fit glmmTMB (with a regularizing ranef prior where needed) to get
## realistic "pretty" (rounded) true parameters, simulate B new datasets
## from those parameters using the real design/covariates, then refit all
## five methods (glmmTMB, joint-phi, PIRLS/digamma-phi, lme4 current
## [moment-phi fix], lme4 2.0-6) to each simulated dataset and record
## status, parameters, derived RE SD/correlation, dispersion, and elapsed
## time.
##
## Unlike the multistart diagnostic (misc/Gamma_GLMM/multistart_*.R), this is
## NOT about starting-point sensitivity for one dataset -- each replicate
## gets one fit per method, from that method's own default/reasonable
## starting values, and the question is parameter recovery and reliability
## (clean/warning/error/singular) across resampled datasets at a range of
## "true" parameter settings.

suppressMessages({
  library(lme4)
  library(reformulas)   # nobars()
})

## ---- joint-phi devfun (generalized from misc/Gamma_GLMM/multistart_jointphi.R) ----
## Builds a devfun that computes -2*logLik directly for a GIVEN phi (no
## internal profiling), over params = c(theta[1:nth], logphi, beta[1:p]).
## Generalized (vs. the multistart-diagnostic version) to arbitrary formula/
## fixed-effect width p, and returns the pieces needed to run bobyqa with
## sensible bounds/starting values for any dataset, not just the one
## two-random-slope toy example.
make_joint_phi_devfun <- function(form, data, family) {
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

    ldL2 <- 2 * determinant(L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL

    nu <- 1 / phi
    ## NOTE: no "+2" here -- unlike Gamma()$aic(), which pads by +2 because
    ## it's designed to be combined with glm's own "+2*rank" elsewhere.
    ## This is a bare -2*logLik term, not an AIC term; the earlier version
    ## of this code copied Gamma()$aic()'s "+2" by mistake, inflating
    ## reported negll by a constant 2 (parameter estimates were unaffected,
    ## since optimizing f+constant has the same optimum, but any numeric
    ## comparison of negll VALUES across methods -- e.g. misc/Gamma_GLMM/multistart_*
    ## -- was off by this constant for the joint-phi arm).
    fit_term <- -2 * sum(weights * dgamma(y, shape = nu, scale = mu * phi, log = TRUE))
    val <- fit_term + sum(u^2) + ldL2
    if (!is.finite(val)) return(1e10)
    val
  }

  ## starting values: lme4's own default theta (diag=1, offdiag=0), a GLM
  ## (fixed effects only, no RE) for beta/dispersion starting values --
  ## the same information lme4 itself uses internally to start PIRLS.
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

## ---- PIRLS with a nested fixed-point phi (digamma-MLE arm) ----
## The outer optimizer (bobyqa) sees only c(theta, beta), exactly like the
## C++ "current" moment-phi fix (glmerLaplace()'s outer fixed-point loop
## over phi is never exposed to the R-level/C++ optimizer either) -- this
## is the R-level twin of that architecture, differing only in which
## estimator re-derives phi from the residual deviance each outer
## iteration: the crude moment plug-in (dev/n, what the C++ fix uses) vs.
## the true Gamma conditional MLE (solving log(nu)-digamma(nu)=dev/(2n)
## for nu=1/phi). This was tried and set aside early in this investigation
## (README #6, "ALSO FALSIFIED") for a single-RE bias-correction question,
## before glmmTMB comparisons or the singular/non-singular multi-RE
## scenarios in this survey existed -- resurrected here as a 5th arm to
## check whether that earlier null result generalizes.
##
## Reuses lme4pureR::pirls(phiType="digamma")'s validated MLE root-finder
## logic (.gammaPhiMLE below is a verbatim copy, since that helper is
## unexported there) but is implemented as its own self-contained closure
## (rather than calling lme4pureR::pirls() directly) so phi can be read
## back out after optimization -- pirls()'s returned function only returns
## the deviance value itself, with no way to recover the internal
## converged phi from outside its closure.
.gammaPhiMLE <- function(dev, n) {
  target <- dev / (2 * n)
  if (!is.finite(target) || target <= 0) return(dev / n)
  f <- function(nu) log(nu) - digamma(nu) - target
  lo <- 1e-8; hi <- 1e8
  if (f(lo) < 0 || f(hi) > 0) return(dev / n)
  1 / uniroot(f, c(lo, hi), tol = 1e-12)$root
}

make_pirls_phi_devfun <- function(form, data, family, phiType = c("moment", "digamma"),
                                   maxPhiIter = 30) {
  phiType <- match.arg(phiType)
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

  last_phi <- NA_real_   ## side-channel: last converged phi, read after optimization

  devfun <- function(thetabeta) {
    theta <- thetabeta[1:nth]
    beta <- thetabeta[(nth + 1):length(thetabeta)]

    Lambdat@x[] <- thfun(theta)
    LtZt <- Lambdat %*% Zt
    offb <- offset + as.vector(X %*% beta)
    eta <- numeric(n); mu <- numeric(n); u <- numeric(q)
    L <- Matrix::Cholesky(Matrix::tcrossprod(LtZt), perm = FALSE, LDL = FALSE, Imult = 1)

    updatemu <- function(uu, phi) {
      eta[] <<- offb + as.vector(crossprod(LtZt, uu))
      mu[] <<- linkinv(eta)
      sum(sqDevResid(y, mu, weights)) / phi + sum(uu^2)
    }
    runPIRLS <- function(phi) {
      u[] <<- numeric(q)
      olducden <- updatemu(u, phi)
      cvgd <- FALSE
      for (i in 1:60) {
        Whalf <- Matrix::Diagonal(x = sqrt(weights / (phi * variance(mu))))
        LtZtMWhalf <- LtZt %*% (Matrix::Diagonal(x = muEta(eta)) %*% Whalf)
        L <<- Matrix::update(L, LtZtMWhalf, 1)
        wtres <- Whalf %*% (y - mu)
        delu <- as.vector(Matrix::solve(L, LtZtMWhalf %*% wtres - u))
        ucden <- updatemu(u + delu, phi)
        if (abs((olducden - ucden) / ucden) < 1e-8) { cvgd <- TRUE; break }
        if (ucden > olducden) {
          for (j in 1:10) {
            ucden <- updatemu(u + (delu <- delu / 2), phi)
            if (ucden <= olducden) break
          }
          if (ucden > olducden) return(FALSE)
        }
        olducden <- ucden
        u[] <<- u + delu
      }
      cvgd
    }

    phi <- 1
    ok <- TRUE
    for (outer in 1:maxPhiIter) {
      ok <- runPIRLS(phi)
      if (!ok) break
      devSum <- sum(sqDevResid(y, mu, weights))
      phiNew <- if (phiType == "moment") devSum / n else .gammaPhiMLE(devSum, n)
      if (abs(phiNew - phi) / phi < 1e-8) { phi <- phiNew; break }
      phi <- phiNew
    }
    if (!ok) return(1e10)
    last_phi <<- phi

    ldL2 <- 2 * determinant(L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL
    nu <- 1 / phi
    ## NOTE: no "+2" here -- unlike Gamma()$aic(), which pads by +2 because
    ## it's designed to be combined with glm's own "+2*rank" elsewhere.
    ## This is a bare -2*logLik term, not an AIC term; the earlier version
    ## of this code copied Gamma()$aic()'s "+2" by mistake, inflating
    ## reported negll by a constant 2 (parameter estimates were unaffected,
    ## since optimizing f+constant has the same optimum, but any numeric
    ## comparison of negll VALUES across methods -- e.g. misc/Gamma_GLMM/multistart_*
    ## -- was off by this constant for the joint-phi arm).
    fit_term <- -2 * sum(weights * dgamma(y, shape = nu, scale = mu * phi, log = TRUE))
    val <- fit_term + sum(u^2) + ldL2
    if (!is.finite(val)) return(1e10)
    val
  }

  glm0 <- glm(nobars(form), data = data, family = family)
  beta_start <- unname(coef(glm0))
  theta_start <- gm$reTrms$theta
  lower_theta <- ifelse(is.finite(gm$reTrms$lower), gm$reTrms$lower, -10)
  upper_theta <- rep(10, nth)

  list(devfun = devfun, nth = nth, p = p, cnms = gm$reTrms$cnms,
       start = c(theta_start, beta_start),
       lower = c(lower_theta, rep(-20, p)),
       upper = c(upper_theta, rep(20, p)),
       beta_names = colnames(X),
       get_phi = function() last_phi)
}

## ---- RE SD/correlation extraction from a raw lme4-native theta vector ----
## Works for any method whose theta uses lme4's raw-Cholesky convention
## (glmer, and the joint-phi/pirls-phi devfuns above, since all are built
## from the same glFormula()/reTrms machinery). Verified (session note)
## that for Gamma GLMMs there is no extra sigma scaling: Sigma_RE =
## Lambda %*% t(Lambda) directly, so sc=1 in mkVarCorr() reproduces
## glmer's own VarCorr() exactly.
##
## Two distinct multi-parameter shapes occur across this survey's
## examples: one grouping factor with >1 correlated term (epil2 complex:
## (Visit|subject), sd has length 2 and a real correlation), or >1
## independent single-term grouping factors (Report4BB: (1|location) +
## (1|fyear), sd has length 2 but there's no correlation between
## different grouping factors to report). Handled by cases below; assumes
## no example in this survey mixes both (a multi-term factor AND more
## than one grouping factor) -- extend further if one is added.
extract_sdcorr <- function(cnms, theta, tol = 1e-5) {
  vc <- lme4:::mkVarCorr(sc = 1, cnms = cnms, theta = theta)
  if (length(vc) == 1) {
    g <- vc[[1]]
    sd <- attr(g, "stddev")
    corr <- attr(g, "correlation")
    list(sd = sd, corr = if (length(sd) > 1) corr[1, 2] else NA_real_,
         singular = det(g) < tol)
  } else {
    sd <- unlist(lapply(vc, function(g) attr(g, "stddev")))
    list(sd = sd, corr = NA_real_,
         singular = any(vapply(vc, det, numeric(1)) < tol))
  }
}

## ---- correlation <-> glmmTMB's unconstrained 2x2 "us" theta transform ----
## Verified numerically (session note): for a 2-RE-term unstructured block,
## glmmTMB's raw correlation theta component satisfies rho = theta/sqrt(1+theta^2).
corr_to_glmmTMB_theta <- function(rho) rho / sqrt(1 - rho^2)

## ---- standardized per-fit result skeleton ----
## Every fitOne_* function below returns this shape so results across the
## four methods can be row-bound directly.
emptyResult <- function(i, status, msg, time_sec = NA_real_,
                         beta = NULL, sd = c(NA_real_, NA_real_), corr = NA_real_,
                         phi = NA_real_, negll = NA_real_, singular = NA) {
  list(i = i, status = status, msg = msg, time_sec = time_sec,
       beta = beta, sd1 = unname(sd[1]), sd2 = unname(sd[2]), corr = corr,
       phi = phi, negll = negll, singular = singular)
}

fit_glmmTMB_one <- function(i, form, dat, family) {
  suppressMessages(library(glmmTMB))
  warn_msgs <- character(0)
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    tryCatch(glmmTMB(form, data = dat, family = family), error = function(e) e),
    warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(fit, "error")) {
    return(emptyResult(i, "error", conditionMessage(fit), time_sec))
  }
  status <- if (length(warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(performance::check_singularity(fit), error = function(e) NA)
  vclist <- VarCorr(fit)$cond
  sc <- if (length(vclist) == 1) {
    vc <- vclist[[1]]; sd <- attr(vc, "stddev"); corr <- attr(vc, "correlation")
    list(sd = sd, corr = if (length(sd) > 1) corr[1, 2] else NA_real_)
  } else {
    list(sd = unlist(lapply(vclist, function(g) attr(g, "stddev"))), corr = NA_real_)
  }
  negll <- tryCatch({
    ll <- as.numeric(logLik(fit))
    if (is.na(ll)) 2 * fit$obj$fn(fit$fit$par) else -2 * ll
  }, error = function(e) NA_real_)
  emptyResult(i, status, paste(warn_msgs, collapse = "; "), time_sec,
              beta = fixef(fit)$cond, sd = sc$sd, corr = sc$corr,
              phi = sigma(fit)^2, negll = negll, singular = singular)
}

fit_glmer_one <- function(i, form, dat, family) {
  warn_msgs <- character(0)
  t0 <- Sys.time()
  fit <- withCallingHandlers(
    tryCatch(glmer(form, data = dat, family = family), error = function(e) e),
    warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(fit, "error")) {
    return(emptyResult(i, "error", conditionMessage(fit), time_sec))
  }
  status <- if (length(warn_msgs) > 0) "warning" else "clean"
  singular <- tryCatch(performance::check_singularity(fit), error = function(e) NA)
  vclist <- VarCorr(fit)
  sc <- if (length(vclist) == 1) {
    vc <- vclist[[1]]; sd <- attr(vc, "stddev"); corr <- attr(vc, "correlation")
    list(sd = sd, corr = if (length(sd) > 1) corr[1, 2] else NA_real_)
  } else {
    list(sd = unlist(lapply(vclist, function(g) attr(g, "stddev"))), corr = NA_real_)
  }
  sd <- sc$sd; corr <- sc$corr
  negll <- tryCatch(-2 * as.numeric(logLik(fit)), error = function(e) NA_real_)
  emptyResult(i, status, paste(warn_msgs, collapse = "; "), time_sec,
              beta = fixef(fit), sd = sd, corr = corr,
              phi = sigma(fit)^2, negll = negll, singular = singular)
}

fit_jointphi_one <- function(i, spec, dat) {
  suppressMessages(library(minqa))
  jp <- make_joint_phi_devfun(spec$formula, dat, spec$family)
  t0 <- Sys.time()
  opt <- tryCatch(
    bobyqa(jp$start, jp$devfun, lower = jp$lower, upper = jp$upper,
           control = list(maxfun = 8000)),
    error = function(e) e
  )
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(opt, "error")) {
    return(emptyResult(i, "error", conditionMessage(opt), time_sec))
  }
  val <- opt$fval
  if (!is.finite(val) || val >= 1e10) {
    return(emptyResult(i, "degenerate", "hit sentinel / non-finite", time_sec))
  }
  status <- if (opt$ierr != 0) "warning" else "clean"
  theta <- opt$par[1:jp$nth]
  phi <- exp(opt$par[jp$nth + 1])
  beta <- opt$par[(jp$nth + 2):length(opt$par)]
  names(beta) <- jp$beta_names
  sc <- extract_sdcorr(jp$cnms, theta)
  emptyResult(i, status, if (opt$ierr != 0) paste("ierr =", opt$ierr) else "",
              time_sec, beta = beta, sd = sc$sd, corr = sc$corr,
              phi = phi, negll = val, singular = sc$singular)
}

fit_pirls_phi_one <- function(i, spec, dat, phiType = "digamma", maxPhiIter = 30) {
  suppressMessages(library(minqa))
  jp <- make_pirls_phi_devfun(spec$formula, dat, spec$family, phiType = phiType,
                               maxPhiIter = maxPhiIter)
  t0 <- Sys.time()
  opt <- tryCatch(
    bobyqa(jp$start, jp$devfun, lower = jp$lower, upper = jp$upper,
           control = list(maxfun = 8000)),
    error = function(e) e
  )
  time_sec <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(opt, "error")) {
    return(emptyResult(i, "error", conditionMessage(opt), time_sec))
  }
  val <- opt$fval
  if (!is.finite(val) || val >= 1e10) {
    return(emptyResult(i, "degenerate", "hit sentinel / non-finite", time_sec))
  }
  status <- if (opt$ierr != 0) "warning" else "clean"
  theta <- opt$par[1:jp$nth]
  beta <- opt$par[(jp$nth + 1):length(opt$par)]
  names(beta) <- jp$beta_names
  sc <- extract_sdcorr(jp$cnms, theta)
  emptyResult(i, status, if (opt$ierr != 0) paste("ierr =", opt$ierr) else "",
              time_sec, beta = beta, sd = sc$sd, corr = sc$corr,
              phi = jp$get_phi(), negll = val, singular = sc$singular)
}

resultsToDF <- function(results, beta_names) {
  status <- vapply(results, `[[`, character(1), "status")
  msg <- vapply(results, `[[`, character(1), "msg")
  time_sec <- vapply(results, `[[`, numeric(1), "time_sec")
  sd1 <- vapply(results, `[[`, numeric(1), "sd1")
  sd2 <- vapply(results, `[[`, numeric(1), "sd2")
  corr <- vapply(results, `[[`, numeric(1), "corr")
  phi <- vapply(results, `[[`, numeric(1), "phi")
  negll <- vapply(results, `[[`, numeric(1), "negll")
  singular <- vapply(results, function(x) isTRUE(x$singular), logical(1))
  beta_mat <- do.call(rbind, lapply(results, function(x) {
    if (is.null(x$beta)) rep(NA_real_, length(beta_names)) else x$beta[beta_names]
  }))
  colnames(beta_mat) <- beta_names
  data.frame(i = seq_along(results), status = status, singular = singular, msg = msg,
             time_sec = time_sec, sd1 = sd1, sd2 = sd2, corr = corr, phi = phi,
             negll = negll, beta_mat, check.names = FALSE)
}
