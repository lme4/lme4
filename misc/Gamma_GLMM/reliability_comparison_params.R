## Follow-up to reliability_comparison.R: for B resimulated datasets
## (identical seed sequence across both lme4 versions), capture the actual
## fitted theta/beta/sigma for every replicate that doesn't error, so we
## can check whether "clean convergence" (no R error/warning) fits are
## actually near the true parameters, or just landed on some other
## non-pathological-looking-but-wrong result.
##
## Parallelized across replicates with parallel::mclapply (each replicate's
## seed is set explicitly and independently inside simulateOne(), so this
## is safe/reproducible regardless of fork scheduling order).

USE_OLD <- TRUE   # TRUE: lme4 2.0-6 (unmodified); FALSE: current dev install
MC_CORES <- 14    # machine has 32 cores; running old+current simultaneously,
                   # 14 each leaves headroom for other work

if (USE_OLD) {
  tmplib <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/lme4_206_lib"
  library(lme4, lib.loc = tmplib)
  stopifnot(packageVersion("lme4") == "2.0-6")
} else {
  library(lme4)
}
library(parallel)
cat("Using lme4", as.character(packageVersion("lme4")),
    if (USE_OLD) "(unmodified, from isolated lib)" else "(current dev install)",
    "-- mc.cores =", MC_CORES, "\n")

simulateOne <- function(seed, n_groups = 20, n_per_group = 20,
                         theta_true = c(0,0,0,1.0,0.5,0.3), beta_true = 2, sigma_true = 2) {
  set.seed(seed)
  n <- n_groups * n_per_group
  dat <- data.frame(
    group1 = rep(1:n_groups, each = n_per_group),
    group2 = rep(1:n_groups, each = n_per_group),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  form2 <- y2 ~ 1 + (1 + x1|group1) + (1 + x2|group2)
  dat$y2 <- simulate(form2[-2], newdata = dat, family = Gamma(link = "log"),
                      newparams = list(theta = theta_true, beta = beta_true, sigma = sigma_true))[[1]]
  list(dat = dat, form2 = form2)
}

theta_true <- c(0,0,0,1.0,0.5,0.3)
beta_true <- 2
sigma_true <- 2

B <- 100
master_seed <- 9000

fitOne <- function(b) {
  sim <- simulateOne(master_seed + b)
  warn_msgs <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      glmer(sim$form2, family = Gamma(link = "log"), data = sim$dat),
      error = function(e) e
    ),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) {
    list(rep = b, status = "error", theta = rep(NA_real_, 6), beta = NA_real_, sigma = NA_real_,
         err_msg = conditionMessage(fit))
  } else {
    status <- if (length(warn_msgs) > 0) "converged with warnings" else "clean convergence"
    list(rep = b, status = status, theta = getME(fit, "theta"), beta = unname(fixef(fit)[1]),
         sigma = sigma(fit), err_msg = NA_character_)
  }
}

t0 <- Sys.time()
results <- mclapply(seq_len(B), fitOne, mc.cores = MC_CORES)
cat("total time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")

status <- vapply(results, function(x) x$status, character(1))
theta_est <- do.call(rbind, lapply(results, function(x) x$theta))
beta_est <- vapply(results, function(x) x$beta, numeric(1))
sigma_est <- vapply(results, function(x) x$sigma, numeric(1))

cat("\n=== true parameters ===\n")
cat("theta:", theta_true, " beta:", beta_true, " sigma:", sigma_true, "\n")

cat("\n=== status counts over", B, "replicates ===\n")
print(table(status))

clean_idx <- status == "clean convergence"
cat("\n=== 'clean convergence' cases (n=", sum(clean_idx), ") -- theta summary ===\n")
print(round(apply(theta_est[clean_idx, , drop=FALSE], 2, summary), 3))
cat("\nbeta summary (true=2):\n"); print(summary(beta_est[clean_idx]))
cat("\nsigma summary (true=2):\n"); print(summary(sigma_est[clean_idx]))

## flag "clean convergence" cases that still landed somewhere clearly wrong:
## group1 diagonal theta (elements 1,3) should be ~0 (true singular); large
## values indicate a pathological but unflagged fit
cat("\n=== 'clean convergence' cases with group1 theta[1] or theta[3] > 0.5",
    "(should be ~0; large values indicate a pathological but unflagged fit) ===\n")
bad_clean <- clean_idx & (theta_est[,1] > 0.5 | theta_est[,3] > 0.5)
cat("count:", sum(bad_clean, na.rm=TRUE), "out of", sum(clean_idx), "clean fits\n")
if (any(bad_clean, na.rm=TRUE)) {
  print(data.frame(rep = which(bad_clean), theta_est[bad_clean, , drop=FALSE],
                    beta = beta_est[bad_clean], sigma = sigma_est[bad_clean]))
}

outfile <- sprintf(
  "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/reliability_params_%s.rds",
  if (USE_OLD) "old" else "current")
saveRDS(list(status = status, theta_est = theta_est, beta_est = beta_est, sigma_est = sigma_est,
             theta_true = theta_true, beta_true = beta_true, sigma_true = sigma_true,
             version = as.character(packageVersion("lme4")), use_old = USE_OLD, B = B,
             err_msgs = vapply(results, function(x) x$err_msg, character(1))),
        outfile)
cat("\nsaved to", outfile, "\n")
