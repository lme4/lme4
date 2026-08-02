## Combine the four multistart runs (lme4 2.0-6, lme4 current/fixed,
## lme4 joint-phi [phi as a first-class optimized parameter], glmmTMB) and
## produce (1) a Raue et al. (2013)-style ECDF ("waterfall") plot of
## achieved -2*logLik across B=200 random starting points each, and (2)
## histograms of the estimated dispersion (phi = sigma^2) across those same
## starts.
##
## NOTE on deviance(fit) vs -2*logLik(fit) for lme4: these are NOT the same
## quantity for these fits (confirmed: e.g. 2195.7 vs 1359.9 on one fit --
## not just a small numerical discrepancy). The lme4 scripts record both
## but use -2*logLik() (field "neg2ll") as the comparable quantity for this
## analysis. glmmTMB's script already used the correct comparable quantity
## (2*fit$obj$fn(fit$fit$par), verified to match -2*logLik(fit) exactly
## when the Hessian is PD) under the field name "deviance" (not rerun here
## since it didn't need correcting). The joint-phi devfun computes this
## quantity directly by construction (aic_term(phi) + ||u||^2 + ldL2 for a
## GIVEN phi, no internal profiling) -- its "neg2ll" field needs no
## correction either.

r_old <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_lme4old.rds")
r_cur <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_lme4current.rds")
r_tmb <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_glmmTMB.rds")
r_joint <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_jointphi.rds")

neg2ll_old <- r_old$neg2ll
neg2ll_cur <- r_cur$neg2ll
neg2ll_tmb <- r_tmb$deviance  ## field name predates the neg2ll rename; same quantity
neg2ll_joint <- r_joint$neg2ll

## a small number of lme4-current fits silently hit the C++ fix's flat 1e10
## sentinel (returned when even the safe phi=1 retry fails inside
## glmerLaplace()) without the outer optimizer raising an R-level error --
## they get an lme4 "warning" status but a nonsensical neg2ll (~1e10).
## Recode these as "degenerate" and exclude from the plots below; they are
## a distinct, already-diagnosed failure mode (misc/Gamma_GLMM/README_Gamma_GLMMs.md
## #9), not part of the multimodality question this script addresses.
sentinel_idx <- which(neg2ll_cur > 1e6)
cat("lme4-current: recoding", length(sentinel_idx), "sentinel-hit fit(s) as 'degenerate' (excluded below)\n")
r_cur$status[sentinel_idx] <- "degenerate"
neg2ll_cur[sentinel_idx] <- NA

cat("\n=== status counts ===\n")
cat("lme4 2.0-6:    "); print(table(r_old$status))
cat("lme4 current:  "); print(table(r_cur$status))
cat("lme4 joint-phi:"); print(table(r_joint$status))
cat("glmmTMB:       "); print(table(r_tmb$status))

cat("\n=== max|deviance(fit) - (-2*logLik(fit))| for lme4 fits ===\n")
cat("old:    ", max(abs(r_old$deviance_fn - r_old$neg2ll), na.rm = TRUE), "\n")
cat("current:", max(abs(r_cur$deviance_fn - r_cur$neg2ll), na.rm = TRUE), "\n")

cat("\n=== distinct -2*logLik values (rounded to nearest 1), with counts ===\n")
cat("lme4 2.0-6:\n"); print(sort(table(round(neg2ll_old[!is.na(neg2ll_old)])), decreasing = TRUE)[1:10])
cat("\nlme4 current:\n"); print(sort(table(round(neg2ll_cur[!is.na(neg2ll_cur)])), decreasing = TRUE)[1:10])
cat("\nlme4 joint-phi:\n"); print(sort(table(round(neg2ll_joint[!is.na(neg2ll_joint)])), decreasing = TRUE)[1:10])
cat("\nglmmTMB:\n"); print(sort(table(round(neg2ll_tmb[!is.na(neg2ll_tmb)])), decreasing = TRUE)[1:10])

## dispersion (phi = sigma^2); sigma() returns sqrt(phi) for both lme4 and glmmTMB
phi_old <- r_old$sigma_est^2
phi_cur <- r_cur$sigma_est^2
phi_cur[sentinel_idx] <- NA
phi_tmb <- r_tmb$sigma_est^2
phi_joint <- r_joint$phi_est
cat("\n=== dispersion (phi = sigma^2) summary, true phi = 4 ===\n")
cat("lme4 2.0-6:    "); print(summary(phi_old))
cat("lme4 current:  "); print(summary(phi_cur))
cat("lme4 joint-phi:"); print(summary(phi_joint))
cat("glmmTMB:       "); print(summary(phi_tmb))

## ECDF plot (corrected: -2*logLik, not deviance(fit)); left panel overlays
## all three lme4 variants (directly comparable -- same devfun family/theta
## units), right panel glmmTMB alone as before
png("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_ecdf.png",
    width = 1000, height = 700, res = 120)
par(mfrow = c(1, 2))

xr <- range(c(neg2ll_old, neg2ll_cur, neg2ll_joint), na.rm = TRUE)
plot(ecdf(neg2ll_old), col = "steelblue", lwd = 2, main = "lme4: old (2.0-6) vs current vs joint-phi",
     xlab = "-2*logLik at convergence", ylab = "ECDF (fraction of starts)",
     xlim = xr, verticals = TRUE, do.points = FALSE)
lines(ecdf(neg2ll_cur), col = "firebrick", lwd = 2, verticals = TRUE, do.points = FALSE)
lines(ecdf(neg2ll_joint), col = "purple", lwd = 2, verticals = TRUE, do.points = FALSE)
legend("bottomright", c("lme4 2.0-6", "lme4 current (fix)", "lme4 joint-phi"),
       col = c("steelblue", "firebrick", "purple"), lwd = 2, bty = "n")

plot(ecdf(neg2ll_tmb), col = "darkgreen", lwd = 2, main = "glmmTMB",
     xlab = "-2*logLik at convergence", ylab = "ECDF (fraction of starts)",
     verticals = TRUE, do.points = FALSE)

dev.off()
cat("\nsaved ECDF plot to multistart_ecdf.png\n")

## dispersion histograms, now four panels
png("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_dispersion_hist.png",
    width = 1500, height = 420, res = 120)
par(mfrow = c(1, 4))
xr_phi <- range(c(phi_old, phi_cur, phi_joint, phi_tmb), na.rm = TRUE)
hist(phi_old, breaks = 30, col = "steelblue", main = "lme4 2.0-6", xlab = "estimated dispersion (phi)", xlim = xr_phi)
abline(v = 4, col = "black", lty = 2, lwd = 2)
hist(phi_cur, breaks = 30, col = "firebrick", main = "lme4 current (fix)", xlab = "estimated dispersion (phi)", xlim = xr_phi)
abline(v = 4, col = "black", lty = 2, lwd = 2)
hist(phi_joint, breaks = 30, col = "purple", main = "lme4 joint-phi", xlab = "estimated dispersion (phi)", xlim = xr_phi)
abline(v = 4, col = "black", lty = 2, lwd = 2)
hist(phi_tmb, breaks = 30, col = "darkgreen", main = "glmmTMB", xlab = "estimated dispersion (phi)", xlim = xr_phi)
abline(v = 4, col = "black", lty = 2, lwd = 2)
dev.off()
cat("saved dispersion histogram to multistart_dispersion_hist.png (dashed line = true phi = 4)\n")
