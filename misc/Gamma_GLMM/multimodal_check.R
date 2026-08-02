## Check whether the "clean but wrong" fits from reliability_comparison_params.R
## (current/fixed lme4, group1 theta silently converging away from its true
## zero) reflect a genuinely multimodal Laplace likelihood surface (two real
## local optima, starting-value dependent) vs. the optimizer just failing to
## find the true-adjacent optimum for some other reason.

library(lme4)

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

## reps where current/fixed lme4 gave a "clean but wrong" fit (from
## reliability_comparison_params_current.log): 14, 15, 26, 32, 48, 55, 61,
## 62, 66, 69, 87, 89, 93. Use a few of these.
master_seed <- 9000
test_reps <- c(14, 26, 48)

for (rep in test_reps) {
  cat("\n\n========== rep", rep, "==========\n")
  sim <- simulateOne(master_seed + rep)

  gm <- glFormula(sim$form2, data = sim$dat, family = Gamma(link = "log"))
  devfun0 <- do.call(mkGlmerDevfun, gm)
  devfun1 <- updateGlmerDevfun(devfun0, gm$reTrms, nAGQ = 1L)

  ## fit 1: default starting values
  fit_default <- glmer(sim$form2, data = sim$dat, family = Gamma(link = "log"))
  theta_default <- getME(fit_default, "theta")
  dev_default <- deviance(fit_default)

  ## fit 2: start AT the true parameters
  fit_true_start <- tryCatch(
    glmer(sim$form2, data = sim$dat, family = Gamma(link = "log"),
          start = list(theta = theta_true, fixef = c(beta_true, rep(0, length(fixef(fit_default)) - 1)))),
    error = function(e) e
  )

  cat("default-start theta:", round(theta_default, 4), " deviance:", dev_default, "\n")
  if (inherits(fit_true_start, "error")) {
    cat("true-param-start: ERROR:", conditionMessage(fit_true_start), "\n")
  } else {
    cat("true-start theta:  ", round(getME(fit_true_start, "theta"), 4), " deviance:", deviance(fit_true_start), "\n")
  }

  ## direct deviance evaluation (no further optimization) at:
  ##  (a) the default-start converged solution
  ##  (b) the true parameters exactly
  beta_default <- unname(fixef(fit_default))
  dev_at_default <- devfun1(c(theta_default, beta_default))
  dev_at_true <- devfun1(c(theta_true, beta_true, rep(0, length(beta_default) - 1)))

  cat("\ndirect devfun evaluation (no optimization):\n")
  cat("  at default-converged point: ", dev_at_default, "\n")
  cat("  at true parameters exactly:", dev_at_true, "\n")
  cat("  (lower deviance = better fit under the Laplace approximation)\n")
  if (dev_at_true < dev_at_default) {
    cat("  ==> TRUE PARAMETERS have LOWER deviance: optimizer failed to find the better/true-adjacent mode\n")
  } else {
    cat("  ==> default-converged point has LOWER (or equal) deviance: genuine second mode, not just an optimizer miss\n")
  }
}
