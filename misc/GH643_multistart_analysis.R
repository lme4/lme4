## Combine the three multistart runs (lme4 2.0-6, lme4 current/fixed,
## glmmTMB) and produce a Raue et al. (2013)-style ECDF ("waterfall") plot
## of achieved deviance across B=200 random starting points each, plus
## numeric summaries.

r_old <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_lme4old.rds")
r_cur <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_lme4current.rds")
r_tmb <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_glmmTMB.rds")

cat("=== status counts ===\n")
cat("lme4 2.0-6:   "); print(table(r_old$status))
cat("lme4 current: "); print(table(r_cur$status))
cat("glmmTMB:      "); print(table(r_tmb$status))

## examine the best (lowest-deviance) mode found by lme4-current, since its
## minimum (2064) is notably below the ~2196 cluster
best_i <- which.min(r_cur$deviance)
cat("\n=== lme4-current best fit found (i=", best_i, ") ===\n")
cat("deviance:", r_cur$deviance[best_i], " status:", r_cur$status[best_i], "\n")
cat("theta:", round(r_cur$theta_est[best_i, ], 4), "\n")
cat("beta:", r_cur$beta_est[best_i], " sigma:", r_cur$sigma_est[best_i], "\n")
cat("start theta:", round(r_cur$start_theta[best_i, ], 4), " start beta:", r_cur$start_beta[best_i], "\n")

## cluster lme4-current's clean-or-warning fits by rounded deviance to see
## how many distinct plateaus/modes are present
cat("\n=== lme4-current: distinct deviance values (rounded to nearest 1), with counts ===\n")
print(sort(table(round(r_cur$deviance[!is.na(r_cur$deviance)])), decreasing = TRUE)[1:10])

cat("\n=== lme4-old: distinct deviance values (rounded to nearest 1), with counts ===\n")
print(sort(table(round(r_old$deviance[!is.na(r_old$deviance)])), decreasing = TRUE)[1:10])

cat("\n=== glmmTMB: distinct deviance values (rounded to nearest 1), with counts ===\n")
print(sort(table(round(r_tmb$deviance[!is.na(r_tmb$deviance)])), decreasing = TRUE)[1:10])

## ECDF plot
png("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_ecdf.png",
    width = 1000, height = 700, res = 120)
par(mfrow = c(1, 2))

## panel 1: lme4 old vs current (directly comparable deviance scale)
xr <- range(c(r_old$deviance, r_cur$deviance), na.rm = TRUE)
plot(ecdf(r_old$deviance), col = "steelblue", lwd = 2, main = "lme4: old (2.0-6) vs current",
     xlab = "deviance (-2logLik) at convergence", ylab = "ECDF (fraction of starts)",
     xlim = xr, verticals = TRUE, do.points = FALSE)
lines(ecdf(r_cur$deviance), col = "firebrick", lwd = 2, verticals = TRUE, do.points = FALSE)
legend("bottomright", c("lme4 2.0-6", "lme4 current (fix)"), col = c("steelblue", "firebrick"), lwd = 2, bty = "n")

## panel 2: glmmTMB on its own scale (not directly comparable to lme4's deviance additive constant)
plot(ecdf(r_tmb$deviance), col = "darkgreen", lwd = 2, main = "glmmTMB",
     xlab = "deviance (2*NLL) at convergence", ylab = "ECDF (fraction of starts)",
     verticals = TRUE, do.points = FALSE)

dev.off()
cat("\nsaved plot to multistart_ecdf.png\n")
