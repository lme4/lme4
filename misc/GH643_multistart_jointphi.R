## Fourth arm for the multistart comparison: phi as a first-class element
## of the parameter vector for the outer nonlinear optimizer (bobyqa),
## rather than profiled via a nested fixed point inside PIRLS (the
## current/fixed lme4 approach) or ignored/fixed at an implicit phi=1
## (old lme4). Adapted from the original prototype misc/GH643_pirls_joint_phi.R
## (which validated this approach for the single-RE bias-correction
## question before the C++ fix existed) onto the same near-singular
## two-random-slope target dataset and B=200 random-start hypercube used
## in misc/GH643_multistart_lme4.R / misc/GH643_multistart_glmmTMB.R, so
## all four methods are directly comparable.
##
## The devfun here directly computes the Laplace deviance
## aic_term(phi) + ||u||^2 + ldL2 for a GIVEN, externally-supplied phi
## (no internal re-estimation at all) -- this literally *is* -2*logLik by
## construction, so there's no deviance()-vs-logLik() ambiguity to worry
## about here (unlike the real lme4 fits; see README #11).

suppressMessages(library(lme4))
library(minqa)   ## bobyqa, lme4's own default optimizer -- fairer comparison
library(parallel)

MC_CORES <- 28   ## matching the other three arms' allocation

## fixed target dataset (seed 9014, matching the other three multistart scripts)
set.seed(9014)
n_groups <- 20; n_per_group <- 20; n <- n_groups * n_per_group
dat <- data.frame(
  group1 = rep(1:n_groups, each = n_per_group),
  group2 = rep(1:n_groups, each = n_per_group),
  x1 = rnorm(n),
  x2 = rnorm(n)
)
form2 <- y2 ~ 1 + (1 + x1|group1) + (1 + x2|group2)
dat$y2 <- simulate(form2[-2], newdata = dat, family = Gamma(link = "log"),
                    newparams = list(theta = c(0,0,0,1,0.5,0.3), beta = 2, sigma = 2))[[1]]

make_devfun_joint <- function(form, data, family) {
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

  ## params = c(theta[1:nth], logphi, beta[1:p])
  function(params) {
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
        if (ucden > olducden) break  ## caller checks cvgd
      }
      olducden <- ucden
      u <- u + delu
    }
    if (!cvgd) return(1e10)  ## same flat-sentinel convention as the C++ fix, for consistency

    ldL2 <- 2 * determinant(L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL

    nu <- 1 / phi
    aic_term <- -2 * sum(weights * dgamma(y, shape = nu, scale = mu * phi, log = TRUE)) + 2
    val <- aic_term + sum(u^2) + ldL2
    if (!is.finite(val)) return(1e10)
    val
  }
}

devfun <- make_devfun_joint(form2, dat, Gamma(link = "log"))
nth <- 6L; p <- 1L

## hypercube: same theta/beta ranges as misc/GH643_multistart_lme4.R, plus
## logphi covering a generous dispersion range around the true phi=4
sampleStart <- function() {
  theta <- numeric(6)
  theta[c(1,3,4,6)] <- runif(4, 0, 2)
  theta[c(2,5)] <- runif(2, -1, 1)
  logphi <- runif(1, log(0.2), log(20))
  beta <- runif(1, 0, 4)
  list(theta = theta, logphi = logphi, beta = beta)
}

B <- 200
set.seed(20260730)  ## same seed as the other three arms
starts <- lapply(seq_len(B), function(i) sampleStart())

## bobyqa bounds: theta diagonal (1,3,4,6) >=0; theta off-diag (2,5), logphi,
## beta unconstrained in practice (wide finite bounds, bobyqa needs finite)
lower <- c(0, -10, 0, 0, -10, 0, log(1e-4), -20)
upper <- c(10, 10, 10, 10, 10, 10, log(1e4), 20)

fitOne <- function(i) {
  s <- starts[[i]]
  start_par <- c(s$theta, s$logphi, s$beta)
  opt <- tryCatch(
    bobyqa(start_par, devfun, lower = lower, upper = upper,
           control = list(maxfun = 8000)),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    return(list(i = i, status = "error", msg = conditionMessage(opt),
                neg2ll = NA_real_, theta = rep(NA_real_, 6), beta = NA_real_, phi = NA_real_))
  }
  val <- opt$fval
  if (!is.finite(val) || val >= 1e10) {
    return(list(i = i, status = "degenerate", msg = "hit sentinel / non-finite",
                neg2ll = NA_real_, theta = rep(NA_real_, 6), beta = NA_real_, phi = NA_real_))
  }
  status <- if (opt$ierr != 0) "warning" else "clean"
  list(i = i, status = status, msg = if (opt$ierr != 0) paste("ierr =", opt$ierr) else NA_character_,
       neg2ll = val, theta = opt$par[1:6], beta = opt$par[8], phi = exp(opt$par[7]))
}

t0 <- Sys.time()
results <- mclapply(seq_len(B), fitOne, mc.cores = MC_CORES)
cat("total time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")

status <- vapply(results, function(x) x$status, character(1))
msg <- vapply(results, function(x) x$msg, character(1))
neg2ll_vec <- vapply(results, function(x) x$neg2ll, numeric(1))
theta_est <- do.call(rbind, lapply(results, function(x) x$theta))
beta_est <- vapply(results, function(x) x$beta, numeric(1))
phi_est <- vapply(results, function(x) x$phi, numeric(1))
sigma_est <- sqrt(phi_est)  ## match the sigma_est convention of the other three arms
start_theta <- do.call(rbind, lapply(starts, function(s) s$theta))
start_beta <- vapply(starts, function(s) s$beta, numeric(1))
start_phi <- vapply(starts, function(s) exp(s$logphi), numeric(1))

cat("\n=== status counts ===\n")
print(table(status))
cat("\nneg2ll summary (excluding NA):\n")
print(summary(neg2ll_vec))
cat("\ndistinct -2*logLik values (rounded to nearest 1), with counts:\n")
print(sort(table(round(neg2ll_vec[!is.na(neg2ll_vec)])), decreasing = TRUE)[1:10])
cat("\nphi summary (true = 4):\n")
print(summary(phi_est))

outfile <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_jointphi.rds"
saveRDS(list(status = status, msg = msg, neg2ll = neg2ll_vec,
             theta_est = theta_est, beta_est = beta_est, sigma_est = sigma_est, phi_est = phi_est,
             start_theta = start_theta, start_beta = start_beta, start_phi = start_phi,
             B = B, seed = 20260730),
        outfile)
cat("\nsaved to", outfile, "\n")
