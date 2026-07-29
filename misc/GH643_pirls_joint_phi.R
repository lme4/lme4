library(lme4)

## Properly-scoped fix: reweight ONLY the deviance/working-weight curvature by
## 1/phi (not the ||u||^2 prior penalty, which is genuinely disp-free since
## the standard GLMM spec has U ~ N(0,I) fixed, NOT N(0,phi*I)). This means
## phi has to be baked into the inner PIRLS mode-finding and working weights,
## not just patched into the final Lm2ll -- so phi becomes a genuine joint
## parameter optimized alongside theta and beta (closer to how glmmTMB does
## full joint ML), rather than a post-hoc dev/n plug-in.

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
  sqDevResid <- fam$dev.resid

  ## params = c(theta, logphi, beta)
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
      ## working weights reweighted by 1/phi (deviance curvature only);
      ## the Imult=1 addition below still adds an UNSCALED +I for the
      ## genuinely disp-free u ~ N(0,I) prior
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
    if (!cvgd) return(NA)

    ldL2 <- 2 * determinant(L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL
    dev_raw <- sum(sqDevResid(y, mu, weights))

    nu <- 1 / phi
    aic_term <- -2 * sum(weights * dgamma(y, shape = nu, scale = mu * phi, log = TRUE)) + 2
    aic_term + sum(u^2) + ldL2
  }
}

fit_devfun_joint <- function(devfun, theta0, logphi0, beta0) {
  start <- c(theta0, logphi0, beta0)
  nth <- length(theta0)
  opt <- optim(start, devfun, method = "Nelder-Mead",
               control = list(maxit = 3000, reltol = 1e-10))
  list(theta = opt$par[seq_len(nth)], phi = exp(opt$par[nth + 1]),
       beta = opt$par[(nth + 2):length(opt$par)], value = opt$value, conv = opt$convergence)
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
  est_theta <- numeric(B); est_phi <- numeric(B)
  for (b in 1:B) {
    b_re <- rnorm(ngrp, 0, sigma_true)
    mu <- exp(mean_ + b_re[dd0$group])
    dd0$y <- rgamma(nrow(dd0), shape = 1/disp, scale = mu*disp)

    df_joint <- make_devfun_joint(form, dd0, Gamma(link="log"))
    r <- tryCatch(fit_devfun_joint(df_joint, theta0 = 1, logphi0 = log(disp), beta0 = 1),
                  error = function(e) NULL)
    est_theta[b] <- if (!is.null(r)) abs(r$theta) else NA
    est_phi[b]   <- if (!is.null(r)) r$phi else NA
  }
  bias <- function(x, truth) 100*(mean(x,na.rm=TRUE)-truth)/truth
  results <- rbind(results, data.frame(
    disp = disp,
    theta_mean = mean(est_theta,na.rm=TRUE), theta_bias = bias(est_theta, sigma_true),
    phi_mean   = mean(est_phi,na.rm=TRUE),   phi_bias   = bias(est_phi, disp),
    nfail = sum(is.na(est_theta))
  ))
  cat(sprintf("disp=%.2f  JOINT: theta bias=%.1f%%   phi bias=%.1f%%  (nfail=%d)\n",
              disp, bias(est_theta, sigma_true), bias(est_phi, disp), sum(is.na(est_theta))))
}

cat("\n=== summary (truth sigma =", sigma_true, ") -- joint (theta,phi,beta) ML ===\n")
print(results)
saveRDS(results, "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/pirls_joint_phi_results.rds")
