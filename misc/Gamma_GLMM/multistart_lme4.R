## Raue et al. (2013) multistart diagnostic: fit the SAME dataset from
## B=200 random starting points drawn from an appropriate hypercube, and
## record -2*logLik() at each (the actually-comparable objective value --
## deviance(fit) is NOT the same quantity for these fits and is recorded
## separately only for reference/diagnostic purposes; see
## misc/Gamma_GLMM/README_Gamma_GLMMs.md #11 update). Plotting the ECDF of -2*logLik
## reveals distinct plateaus if there are multiple genuinely separated
## local optima (multimodality), vs. a single dominant plateau if the
## optimizer reliably finds one (likely global) optimum regardless of
## starting point.
##
## Same target dataset used throughout this investigation: rep=14 from the
## B=100 reliability comparison (seed = master_seed(9000) + 14 = 9014),
## one of the "clean but wrong" cases where the current/fixed lme4 finds a
## non-singular mode for group1 even starting from the true parameters.

USE_OLD <- TRUE   # TRUE: lme4 2.0-6 (unmodified); FALSE: current dev install
MC_CORES <- 10    # 32 cores total; running old lme4 + current lme4 + glmmTMB
                   # simultaneously, ~10 each leaves headroom

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

## fixed target dataset (seed 9014, matching misc/Gamma_GLMM/multimodal_check.R rep=14)
set.seed(9014)
n_groups <- 20; n_per_group <- 20; n <- n_groups * n_per_group
dat <- data.frame(
  group1 = rep(1:n_groups, each = n_per_group),
  group2 = rep(1:n_groups, each = n_per_group),
  x1 = rnorm(n),
  x2 = rnorm(n)
)
form2 <- y2 ~ 1 + (1 + x1|group1) + (1 + x2|group2)
dat$y2 <- simulate(form2[-2], newdata = dat, family = Gamma(link = "log"),
                    newparams = list(theta = c(0,0,0,1,0.5,0.3), beta = 2, sigma = 2))[[1]]

## hypercube for random starting points, in lme4's own theta parameterization
## (raw relative-Cholesky-factor entries; diagonal >=0, off-diagonal free):
##   theta[1], theta[3] (group1 diag), theta[4], theta[6] (group2 diag): U(0,2)
##   theta[2] (group1 offdiag), theta[5] (group2 offdiag): U(-1,1)
##   beta (intercept, only fixed effect): U(0,4)
sampleStart <- function() {
  theta <- numeric(6)
  theta[c(1,3,4,6)] <- runif(4, 0, 2)
  theta[c(2,5)] <- runif(2, -1, 1)
  beta <- runif(1, 0, 4)
  list(theta = theta, beta = beta)
}

B <- 200
set.seed(20260730)
starts <- lapply(seq_len(B), function(i) sampleStart())

fitOne <- function(i) {
  s <- starts[[i]]
  warn_msgs <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      glmer(form2, data = dat, family = Gamma(link = "log"),
            start = list(theta = s$theta, fixef = s$beta)),
      error = function(e) e
    ),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) {
    return(list(i = i, status = "error", msg = conditionMessage(fit),
                neg2ll = NA_real_, deviance_fn = NA_real_,
                theta = rep(NA_real_, 6), beta = NA_real_, sigma = NA_real_))
  }
  status <- if (length(warn_msgs) > 0) "warning" else "clean"
  ll <- as.numeric(logLik(fit))
  list(i = i, status = status, msg = paste(warn_msgs, collapse = "; "),
       neg2ll = -2 * ll, deviance_fn = deviance(fit),
       theta = getME(fit, "theta"), beta = unname(fixef(fit)[1]), sigma = sigma(fit))
}

t0 <- Sys.time()
results <- mclapply(seq_len(B), fitOne, mc.cores = MC_CORES)
cat("total time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")

status <- vapply(results, function(x) x$status, character(1))
msg <- vapply(results, function(x) x$msg, character(1))
neg2ll_vec <- vapply(results, function(x) x$neg2ll, numeric(1))
deviance_fn_vec <- vapply(results, function(x) x$deviance_fn, numeric(1))
theta_est <- do.call(rbind, lapply(results, function(x) x$theta))
beta_est <- vapply(results, function(x) x$beta, numeric(1))
sigma_est <- vapply(results, function(x) x$sigma, numeric(1))
start_theta <- do.call(rbind, lapply(starts, function(s) s$theta))
start_beta <- vapply(starts, function(s) s$beta, numeric(1))

cat("\n=== status counts ===\n")
print(table(status))
cat("\n-2*logLik() summary (excluding NA) -- the comparable quantity:\n")
print(summary(neg2ll_vec))
cat("\ndeviance(fit) summary (excluding NA) -- NOT the same quantity, reference only:\n")
print(summary(deviance_fn_vec))
cat("\nmax|deviance(fit) - (-2*logLik(fit))| across all fits:",
    max(abs(deviance_fn_vec - neg2ll_vec), na.rm = TRUE),
    " (if this varies across fits rather than being ~constant, the two are not just a fixed additive offset apart)\n")

outfile <- sprintf(
  "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_%s.rds",
  if (USE_OLD) "lme4old" else "lme4current")
saveRDS(list(status = status, msg = msg, neg2ll = neg2ll_vec, deviance_fn = deviance_fn_vec,
             theta_est = theta_est, beta_est = beta_est, sigma_est = sigma_est,
             start_theta = start_theta, start_beta = start_beta,
             version = as.character(packageVersion("lme4")), use_old = USE_OLD, B = B,
             seed = 20260730),
        outfile)
cat("\nsaved to", outfile, "\n")
