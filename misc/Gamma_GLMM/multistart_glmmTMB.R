## Raue et al. (2013) multistart diagnostic for glmmTMB, matching
## multistart_lme4.R: same target dataset (seed 9014), B=200 random starts.
## Deviance is retrieved as 2*fit$obj$fn(fit$fit$par) rather than
## -2*logLik(fit), since logLik() returns NA when the Hessian is
## non-positive-definite (expected for singular fits) -- per
## https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html
## fit$obj$fn() still gives the objective (NLL) at the fitted parameters
## regardless of Hessian status. Verified 2*fit$obj$fn(fit$fit$par) matches
## -2*logLik(fit) exactly on a normally-converged fit before trusting this
## at scale.

MC_CORES <- 10   # 32 cores total; running old lme4 + current lme4 + glmmTMB
                  # simultaneously, ~10 each leaves headroom

suppressMessages(library(lme4))  ## only for simulate(); not used for fitting
library(glmmTMB)
library(parallel)
cat("Using glmmTMB", as.character(packageVersion("glmmTMB")), "-- mc.cores =", MC_CORES, "\n")

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

## hypercube for random starting points, in glmmTMB's own (native, default
## zero-initialized) theta parameterization -- log-SD/unconstrained-corr
## mix, not directly comparable to lme4's raw-Cholesky theta, but a broad
## symmetric range around the zero default is the natural analog:
##   theta[1:6] ~ U(-2,2); beta (intercept) ~ U(0,4) (same range as lme4,
##   since beta IS directly comparable -- both are the log-link intercept)
sampleStart <- function() {
  list(theta = runif(6, -2, 2), beta = runif(1, 0, 4))
}

B <- 200
set.seed(20260730)
starts <- lapply(seq_len(B), function(i) sampleStart())

fitOne <- function(i) {
  s <- starts[[i]]
  warn_msgs <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      glmmTMB(form2, data = dat, family = Gamma(link = "log"),
              start = list(theta = s$theta, beta = s$beta)),
      error = function(e) e
    ),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) {
    return(list(i = i, status = "error", msg = conditionMessage(fit),
                deviance = NA_real_, logLik = NA_real_,
                theta_native = rep(NA_real_, 6), beta = NA_real_, sigma = NA_real_))
  }
  dev <- tryCatch(2 * fit$obj$fn(fit$fit$par), error = function(e) NA_real_)
  ll <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)  ## NA if Hessian non-PD; deviance/dev above is the robust one
  status <- if (length(warn_msgs) > 0) "warning" else "clean"
  fitpar <- fit$fit$par
  list(i = i, status = status, msg = paste(warn_msgs, collapse = "; "),
       deviance = dev, logLik = ll,
       theta_native = unname(fitpar[names(fitpar) == "theta"]),
       beta = unname(fitpar[names(fitpar) == "beta"][1]),
       sigma = tryCatch(sigma(fit), error = function(e) NA_real_))
}

t0 <- Sys.time()
results <- mclapply(seq_len(B), fitOne, mc.cores = MC_CORES)
cat("total time:", as.numeric(Sys.time() - t0, units = "mins"), "min\n")

status <- vapply(results, function(x) x$status, character(1))
msg <- vapply(results, function(x) x$msg, character(1))
deviance_vec <- vapply(results, function(x) x$deviance, numeric(1))
logLik_vec <- vapply(results, function(x) x$logLik, numeric(1))
theta_est <- do.call(rbind, lapply(results, function(x) x$theta_native))
beta_est <- vapply(results, function(x) x$beta, numeric(1))
sigma_est <- vapply(results, function(x) x$sigma, numeric(1))
start_theta <- do.call(rbind, lapply(starts, function(s) s$theta))
start_beta <- vapply(starts, function(s) s$beta, numeric(1))

cat("\n=== status counts ===\n")
print(table(status))
cat("\ndeviance summary (excluding NA):\n")
print(summary(deviance_vec))

outfile <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/multistart_glmmTMB.rds"
saveRDS(list(status = status, msg = msg, deviance = deviance_vec, logLik = logLik_vec,
             theta_est = theta_est, beta_est = beta_est, sigma_est = sigma_est,
             start_theta = start_theta, start_beta = start_beta,
             version = as.character(packageVersion("glmmTMB")), B = B, seed = 20260730),
        outfile)
cat("\nsaved to", outfile, "\n")
