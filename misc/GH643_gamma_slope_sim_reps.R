library(lme4)

n <- 1e4
f <- factor(rep(1:20, each = n / 20))

beta_true  <- c(1, 0.5)
sd_int_true  <- 0.3
sd_slope_true <- 0.2
corr_true <- 0
shape_true <- 0.5
sigma_true <- 1 / sqrt(shape_true)
theta_true <- c(sd_int_true, corr_true * sd_slope_true, sd_slope_true)

form <- y ~ 1 + x + (1 + x | f)
B <- 40

res <- data.frame(rep=1:B, b0=NA, b1=NA, sd_int=NA, sd_slope=NA, corr=NA, shape_fit=NA)
nfail <- 0

set.seed(2027)
for (i in 1:B) {
  x <- rnorm(n)
  dd <- data.frame(f = f, x = x)
  dd$y <- simulate(~ 1 + x + (1 + x | f), newdata = dd, family = Gamma(link = "log"),
                    newparams = list(beta = beta_true, theta = theta_true, sigma = sigma_true))[[1]]
  m <- tryCatch(suppressWarnings(glmer(form, data = dd, family = Gamma(link = "log"))),
                error = function(e) NULL)
  if (is.null(m)) { nfail <- nfail + 1; next }
  vc <- as.data.frame(VarCorr(m))
  res$b0[i] <- fixef(m)[1]
  res$b1[i] <- fixef(m)[2]
  res$sd_int[i] <- vc$sdcor[vc$var1=="(Intercept)" & is.na(vc$var2)]
  res$sd_slope[i] <- vc$sdcor[vc$var1=="x" & is.na(vc$var2)]
  res$corr[i] <- vc$sdcor[!is.na(vc$var2)]
  res$shape_fit[i] <- 1/sigma(m)^2
  cat(sprintf("rep %d: sd_int=%.3f sd_slope=%.3f corr=%.3f b0=%.3f b1=%.3f shape=%.3f\n",
              i, res$sd_int[i], res$sd_slope[i], res$corr[i], res$b0[i], res$b1[i], res$shape_fit[i]))
}

cat("\nfailed fits:", nfail, "\n\n")
bias <- function(est, truth) 100*(mean(est,na.rm=TRUE)-truth)/truth
cat("=== summary over", B, "reps (shape_true =", shape_true, ", i.e. disp =", 1/shape_true, ") ===\n")
cat(sprintf("b0:       true=%.3f  mean_est=%.3f  bias=%.1f%%\n", beta_true[1], mean(res$b0,na.rm=T), bias(res$b0,beta_true[1])))
cat(sprintf("b1:       true=%.3f  mean_est=%.3f  bias=%.1f%%\n", beta_true[2], mean(res$b1,na.rm=T), bias(res$b1,beta_true[2])))
cat(sprintf("sd_int:   true=%.3f  mean_est=%.3f  bias=%.1f%%\n", sd_int_true, mean(res$sd_int,na.rm=T), bias(res$sd_int,sd_int_true)))
cat(sprintf("sd_slope: true=%.3f  mean_est=%.3f  bias=%.1f%%\n", sd_slope_true, mean(res$sd_slope,na.rm=T), bias(res$sd_slope,sd_slope_true)))
cat(sprintf("corr:     true=%.3f  mean_est=%.3f\n", corr_true, mean(res$corr,na.rm=T)))
cat(sprintf("shape:    true=%.3f  mean_est=%.3f  bias=%.1f%%\n", shape_true, mean(res$shape_fit,na.rm=T), bias(res$shape_fit,shape_true)))

saveRDS(res, "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/gamma_slope_sim_reps.rds")
