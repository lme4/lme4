library(lme4)

## Third variant: replace the crude phi_hat = dev/n moment-based plug-in
## with the true (Gamma-specific) conditional MLE of the dispersion, found
## by solving  log(nu) - digamma(nu) = dev/(2n)  for the shape nu = 1/phi
## (standard Gamma-GLM dispersion MLE; e.g. MASS::gamma.shape). This equation
## is monotonic in nu (LHS: +inf at nu->0, ->0 as nu->inf), so uniroot finds
## a unique solution.

phi_mle_gamma <- function(dev, W) {
  target <- dev / (2 * W)
  if (!is.finite(target) || target <= 0) return(dev / W)
  f <- function(nu) log(nu) - digamma(nu) - target
  lo <- 1e-8; hi <- 1e8
  if (f(lo) < 0 || f(hi) > 0) return(dev / W)  ## fallback to moment estimator
  1 / uniroot(f, c(lo, hi), tol = 1e-12)$root
}

aic_gamma_phi <- function(y, mu, wt, phi) {
  nu <- 1 / phi
  -2 * sum(wt * dgamma(y, shape = nu, scale = mu * phi, log = TRUE)) + 2
}

make_devfun <- function(form, data, family, variant = c("orig", "momentfix", "mle")) {
  variant <- match.arg(variant)
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
  betaind <- -seq_len(nth)

  fam <- if (is.function(family)) family() else family
  linkinv <- fam$linkinv; variance <- fam$variance; muEta <- fam$mu.eta
  aicfun <- fam$aic; sqDevResid <- fam$dev.resid

  function(thetabeta) {
    Lambdat@x[] <- thfun(thetabeta[-betaind])
    LtZt <- Lambdat %*% Zt
    beta <- thetabeta[betaind]
    offb <- offset + as.vector(X %*% beta)
    eta <- numeric(n); mu <- numeric(n)
    updatemu <- function(uu) {
      eta[] <<- offb + as.vector(crossprod(LtZt, uu))
      mu[] <<- linkinv(eta)
      sum(sqDevResid(y, mu, weights)) + sum(uu^2)
    }
    u <- numeric(q)
    olducden <- updatemu(u)
    L <- Matrix::Cholesky(Matrix::tcrossprod(LtZt), perm = FALSE, LDL = FALSE, Imult = 1)
    cvgd <- FALSE
    for (i in 1:60) {
      Whalf <- Matrix::Diagonal(x = sqrt(weights / variance(mu)))
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
        if (ucden > olducden) stop("step-halving failed")
      }
      olducden <- ucden
      u <- u + delu
    }
    if (!cvgd) return(NA)

    ldL2 <- 2 * determinant(L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL
    dev_raw <- sum(sqDevResid(y, mu, weights))

    if (variant == "orig") {
      aicfun(y, rep.int(1, n), mu, weights, dev_raw) + sum(u^2) + ldL2
    } else if (variant == "momentfix") {
      phi_hat <- dev_raw / n
      aicfun(y, rep.int(1, n), mu, weights, dev_raw) + sum(u^2) / phi_hat + ldL2
    } else {  ## mle
      phi_hat <- phi_mle_gamma(dev_raw, n)
      aic_gamma_phi(y, mu, weights, phi_hat) + sum(u^2) / phi_hat + ldL2
    }
  }
}

fit_devfun <- function(devfun, theta0, beta0) {
  start <- c(theta0, beta0)
  nth <- length(theta0)
  opt <- optim(start, devfun, method = "Nelder-Mead",
               control = list(maxit = 2000, reltol = 1e-10))
  list(theta = opt$par[seq_len(nth)], beta = opt$par[-seq_len(nth)], value = opt$value, conv = opt$convergence)
}

## quick sanity check: does phi_mle_gamma recover known dispersion accurately?
set.seed(1)
yy <- rgamma(2000, shape = 1/0.05, scale = 2*0.05)
devtest <- sum(Gamma()$dev.resid(yy, rep(mean(yy), length(yy)), rep(1,length(yy))))
cat("sanity: true disp=0.05, moment est=", devtest/length(yy),
    " mle est=", phi_mle_gamma(devtest, length(yy)), "\n\n")

## --- single-RE test matching earlier disp_scan.R / pirls_dispfix.R conditions ---
set.seed(2026)
ngrp <- 30; nrep <- 20
sigma_true <- 0.3
mean_ <- 1
dispvals <- c(0.05, 0.2, 1)
B <- 25

dd0 <- expand.grid(group = factor(1:ngrp), rep = 1:nrep)
form <- y ~ 1 + (1 | group)

results <- data.frame()
for (disp in dispvals) {
  est_orig <- numeric(B); est_mfix <- numeric(B); est_mle <- numeric(B)
  for (b in 1:B) {
    b_re <- rnorm(ngrp, 0, sigma_true)
    mu <- exp(mean_ + b_re[dd0$group])
    dd0$y <- rgamma(nrow(dd0), shape = 1/disp, scale = mu*disp)

    df_orig <- make_devfun(form, dd0, Gamma(link="log"), "orig")
    df_mfix <- make_devfun(form, dd0, Gamma(link="log"), "momentfix")
    df_mle  <- make_devfun(form, dd0, Gamma(link="log"), "mle")

    r0 <- tryCatch(fit_devfun(df_orig, theta0 = 1, beta0 = c(1)), error = function(e) NULL)
    r1 <- tryCatch(fit_devfun(df_mfix, theta0 = 1, beta0 = c(1)), error = function(e) NULL)
    r2 <- tryCatch(fit_devfun(df_mle,  theta0 = 1, beta0 = c(1)), error = function(e) NULL)

    est_orig[b] <- if (!is.null(r0)) abs(r0$theta) else NA
    est_mfix[b] <- if (!is.null(r1)) abs(r1$theta) else NA
    est_mle[b]  <- if (!is.null(r2)) abs(r2$theta) else NA
  }
  bias <- function(x) 100*(mean(x,na.rm=TRUE)-sigma_true)/sigma_true
  results <- rbind(results, data.frame(
    disp = disp,
    orig_mean = mean(est_orig,na.rm=TRUE),  orig_bias = bias(est_orig),
    mfix_mean = mean(est_mfix,na.rm=TRUE),  mfix_bias = bias(est_mfix),
    mle_mean  = mean(est_mle,na.rm=TRUE),   mle_bias  = bias(est_mle),
    nfail_orig = sum(is.na(est_orig)), nfail_mfix = sum(is.na(est_mfix)), nfail_mle = sum(is.na(est_mle))
  ))
  cat(sprintf("disp=%.2f  orig bias=%.1f%%   momentfix bias=%.1f%%   MLE-fix bias=%.1f%%\n",
              disp, bias(est_orig), bias(est_mfix), bias(est_mle)))
}

cat("\n=== summary (truth sigma =", sigma_true, ") ===\n")
print(results)
saveRDS(results, "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/pirls_dispfix_mle_results.rds")
