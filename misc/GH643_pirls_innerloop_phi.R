library(lme4)

## Less invasive version of the fix: phi is profiled out via a fixed-point
## loop NESTED INSIDE the PIRLS/devfun evaluation, exactly like beta already
## is for nAGQ=0. The OUTER optimizer only ever sees (theta, beta) -- same
## interface glmer already has -- no new free parameter, no change to the
## outer optim()/bobyqa call. This directly tests whether the architecture
## concern (promoting phi to a genuine outer ML parameter) can be avoided.

phi_mle_gamma <- function(dev, W) {
  target <- dev / (2 * W)
  if (!is.finite(target) || target <= 0) return(dev / W)
  f <- function(nu) log(nu) - digamma(nu) - target
  lo <- 1e-8; hi <- 1e8
  if (f(lo) < 0 || f(hi) > 0) return(dev / W)
  1 / uniroot(f, c(lo, hi), tol = 1e-12)$root
}

make_devfun_inner <- function(form, data, family) {
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
  sqDevResid <- fam$dev.resid

  pirls_given_phi <- function(beta, phi) {
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
        if (ucden > olducden) stop("step-halving failed")
      }
      olducden <- ucden
      u <- u + delu
    }
    list(u = u, mu = mu, L = L, cvgd = cvgd)
  }

  LtZt <- NULL  ## set per-call once theta is known

  function(thetabeta) {
    theta <- thetabeta[-betaind]
    beta <- thetabeta[betaind]
    Lambdat@x[] <<- thfun(theta)
    LtZt <<- Lambdat %*% Zt

    ## nested fixed point over phi: alternate PIRLS-given-phi with
    ## re-estimating phi (Gamma MLE) from the resulting fit, to convergence.
    ## outer optimizer never sees phi -- interface is unchanged (theta,beta).
    phi <- 0.5  ## generic starting guess
    for (outer_it in 1:30) {
      fit <- pirls_given_phi(beta, phi)
      if (!fit$cvgd) return(NA)
      dev_raw <- sum(sqDevResid(y, fit$mu, weights))
      phi_new <- phi_mle_gamma(dev_raw, n)
      if (abs(phi_new - phi) / phi < 1e-8) { phi <- phi_new; break }
      phi <- phi_new
    }

    ldL2 <- 2 * determinant(fit$L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL
    nu <- 1 / phi
    aic_term <- -2 * sum(weights * dgamma(y, shape = nu, scale = fit$mu * phi, log = TRUE)) + 2
    aic_term + sum(fit$u^2) + ldL2
  }
}

fit_devfun <- function(devfun, theta0, beta0) {
  start <- c(theta0, beta0)
  nth <- length(theta0)
  opt <- optim(start, devfun, method = "Nelder-Mead",
               control = list(maxit = 2000, reltol = 1e-10))
  list(theta = opt$par[seq_len(nth)], beta = opt$par[-seq_len(nth)], value = opt$value, conv = opt$convergence)
}

## --- same single-RE test conditions as before ---
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
  est_theta <- numeric(B)
  for (b in 1:B) {
    b_re <- rnorm(ngrp, 0, sigma_true)
    mu <- exp(mean_ + b_re[dd0$group])
    dd0$y <- rgamma(nrow(dd0), shape = 1/disp, scale = mu*disp)

    df_inner <- make_devfun_inner(form, dd0, Gamma(link="log"))
    r <- tryCatch(fit_devfun(df_inner, theta0 = 1, beta0 = 1), error = function(e) NULL)
    est_theta[b] <- if (!is.null(r)) abs(r$theta) else NA
  }
  bias <- function(x) 100*(mean(x,na.rm=TRUE)-sigma_true)/sigma_true
  results <- rbind(results, data.frame(disp = disp, theta_mean = mean(est_theta,na.rm=TRUE),
                                        theta_bias = bias(est_theta), nfail = sum(is.na(est_theta))))
  cat(sprintf("disp=%.2f  INNER-LOOP-PHI (outer sees only theta,beta): theta bias=%.1f%%  (nfail=%d)\n",
              disp, bias(est_theta), sum(is.na(est_theta))))
}

cat("\n=== summary (truth sigma =", sigma_true, ") -- phi profiled inside PIRLS, outer interface unchanged ===\n")
print(results)
saveRDS(results, "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/pirls_innerloop_phi_results.rds")
