## Pairs plots of theta (6 components, native parameterization per engine),
## beta (intercept), and sigma (dispersion sqrt) across the B=200 multistart
## fits, one plot per method (lme4 2.0-6, lme4 current, glmmTMB).

r_old <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_lme4old.rds")
r_cur <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_lme4current.rds")
r_tmb <- readRDS("/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_glmmTMB.rds")

## recode the 2 lme4-current sentinel-hit fits as missing (see #11 in README)
sentinel_idx <- which(r_cur$neg2ll > 1e6)
r_cur$theta_est[sentinel_idx, ] <- NA
r_cur$beta_est[sentinel_idx] <- NA
r_cur$sigma_est[sentinel_idx] <- NA

makeDF <- function(r) {
  df <- data.frame(r$theta_est)
  names(df) <- paste0("theta", 1:6)
  df$beta <- r$beta_est
  df$sigma <- r$sigma_est
  df[complete.cases(df), ]
}

df_old <- makeDF(r_old)
df_cur <- makeDF(r_cur)
df_tmb <- makeDF(r_tmb)

cat("n complete cases: old =", nrow(df_old), " current =", nrow(df_cur), " glmmTMB =", nrow(df_tmb), "\n")

outdir <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb"

png(file.path(outdir, "multistart_pairs_lme4old.png"), width = 1400, height = 1400, res = 130)
pairs(df_old, gap = 0, main = sprintf("lme4 2.0-6: theta (native), beta, sigma (n=%d)", nrow(df_old)))
dev.off()

png(file.path(outdir, "multistart_pairs_lme4current.png"), width = 1400, height = 1400, res = 130)
pairs(df_cur, gap = 0, main = sprintf("lme4 current (fix): theta (native), beta, sigma (n=%d)", nrow(df_cur)))
dev.off()

png(file.path(outdir, "multistart_pairs_glmmTMB.png"), width = 1400, height = 1400, res = 130)
pairs(df_tmb, gap = 0, main = sprintf("glmmTMB: theta (native), beta, sigma (n=%d)", nrow(df_tmb)))
dev.off()

cat("saved multistart_pairs_lme4old.png, multistart_pairs_lme4current.png, multistart_pairs_glmmTMB.png\n")
