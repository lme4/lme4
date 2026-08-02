library(lme4)

## exact design from yesterday's disp_scan.R (single RE, 30 groups x 20 reps
## = 600 obs, true sigma=0.3), now specifically at shape=0.5 (disp=2), to
## isolate whether the random-slope/correlation structure was masking the
## bias in gamma_slope_sim_reps.R

set.seed(2026)
ngrp <- 30; nrep <- 20
sigma_true <- 0.3
mean_ <- 1
disp <- 2          # shape = 0.5
B <- 100

dd0 <- expand.grid(group = factor(1:ngrp), rep = 1:nrep)

est <- numeric(B)
nfail <- 0
for (b in 1:B) {
  b_re <- rnorm(ngrp, 0, sigma_true)
  mu <- exp(mean_ + b_re[dd0$group])
  dd0$y <- rgamma(nrow(dd0), shape = 1/disp, scale = mu*disp)
  m <- tryCatch(suppressWarnings(glmer(y ~ (1|group), data=dd0, family=Gamma(link="log"))),
                error=function(e) NULL)
  if (is.null(m)) { nfail <- nfail + 1; est[b] <- NA; next }
  est[b] <- sqrt(as.numeric(VarCorr(m)$group))
}

cat("failed fits:", nfail, "\n")
m_est <- mean(est, na.rm=TRUE)
cat(sprintf("shape=0.5 (disp=2), single RE, 30 groups x 20 obs: true sigma=%.2f  mean est=%.4f  pct bias=%.1f%%\n",
            sigma_true, m_est, 100*(m_est-sigma_true)/sigma_true))
cat("quantiles of estimates:\n")
print(quantile(est, c(0,0.1,0.25,0.5,0.75,0.9,1), na.rm=TRUE))
