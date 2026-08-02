library(lme4)

## test prediction: glmer's Laplace RE-variance bias should grow as true
## dispersion moves away from 1 (shape moves away from 1), since the
## ldL2/precision term in the Laplace objective is disp-free while the
## aic()/deviance term is not -- an inconsistency that should vanish at disp=1

set.seed(2026)
ngrp <- 30; nrep <- 20
sigma_true <- 0.3
mean_ <- 1
dispvals <- c(0.05, 0.2, 1, 3, 8)
B <- 40

dd0 <- expand.grid(group = factor(1:ngrp), rep = 1:nrep)

results <- data.frame()
for (disp in dispvals) {
  ests <- numeric(B)
  nfail <- 0
  for (b in 1:B) {
    b_re <- rnorm(ngrp, 0, sigma_true)
    mu <- exp(mean_ + b_re[dd0$group])
    dd0$y <- rgamma(nrow(dd0), shape = 1/disp, scale = mu*disp)
    m <- tryCatch(suppressWarnings(glmer(y ~ (1|group), data = dd0, family = Gamma(link="log"))),
                  error = function(e) NULL)
    if (is.null(m)) { nfail <- nfail + 1; ests[b] <- NA; next }
    ests[b] <- sqrt(as.numeric(VarCorr(m)$group))
  }
  est_mean <- mean(ests, na.rm=TRUE)
  results <- rbind(results, data.frame(disp=disp, shape=1/disp,
                                        est_sd=est_mean,
                                        pct_bias=100*(est_mean-sigma_true)/sigma_true,
                                        nfail=nfail))
  cat(sprintf("disp=%.2f (shape=%.2f): mean est sd(RE) = %.4f  pct_bias = %.1f%%  (nfail=%d)\n",
              disp, 1/disp, est_mean, 100*(est_mean-sigma_true)/sigma_true, nfail))
}

cat("\n=== summary (truth sigma =", sigma_true, ") ===\n")
print(results)
saveRDS(results, "disp_scan_results.rds")
