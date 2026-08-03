## Reprex: does plain PIRLS already fail for most (theta,beta) in a
## neighborhood of the near-singular truth used in lme4's
## test-isSingular.R "checking singular fit for merMod" scenario,
## independent of any optimizer trajectory? Evaluates the Laplace deviance
## directly at many perturbed (theta,beta) points (rather than running a
## full glmer() optimization), to isolate PIRLS's own robustness from how
## any particular optimizer happens to walk through parameter space.
##
## Set USE_OLD to compare unmodified CRAN-era lme4 (2.0-6, installed via
## remotes::install_version into an isolated library so it doesn't clobber
## the current dev install) against the current (dev, Gamma_GLMM branch)
## lme4 on the *same* perturbation grid.

USE_OLD <- FALSE   # TRUE: lme4 2.0-6 (unmodified); FALSE: current dev install

if (USE_OLD) {
  tmplib <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/lme4_206_lib"
  library(lme4, lib.loc = tmplib)
  stopifnot(packageVersion("lme4") == "2.0-6")
} else {
  library(lme4)
}
cat("Using lme4", as.character(packageVersion("lme4")),
    if (USE_OLD) paste("(unmodified, from isolated lib)") else "(current dev install)", "\n")

## reconstruct the exact dataset from tests/testthat/test-isSingular.R
## ("checking singular fit for merMod")
set.seed(101)
n_groups <- 20; n_per_group <- 20; n <- n_groups * n_per_group
dat <- data.frame(
  group1 = rep(1:n_groups, each = n_per_group),
  group2 = rep(1:n_groups, each = n_per_group),
  x1 = rnorm(n),
  x2 = rnorm(n)
)
form <- y ~ 1 + (1 + x1|group1) + (1 + x2|group2)
dat$y <- simulate(form[-2], newdata = dat, family = gaussian,
                   newparams = list(beta = c(-2), theta = c(2,3,4,2,3,4), sigma = 2))[[1]]
mod <- lmer(form, data = dat)

form2 <- y2 ~ 1 + (1 + x1|group1) + (1 + x2|group2)
dat$y2 <- simulate(form[-2], newdata = dat, family = gaussian,
                    newparams = list(beta = c(-2), theta = c(0,0,0,2.5,1.5,0.8), sigma = 2))[[1]]
mod2 <- suppressMessages(suppressWarnings(lmer(form2, data = dat)))

form <- y ~ 1 + (1 + x1|group1) + (1 + x2|group2)
dat$y <- simulate(form[-2], newdata = dat, family = poisson(link = "log"),
                   newparams = list(theta = c(0,0,0,1.0,0.5,0.3), beta = 2))[[1]]

form2 <- y2 ~ 1 + (1 + x1|group1) + (1 + x2|group2)
dat$y2 <- simulate(form2[-2], newdata = dat, family = Gamma(link = "log"),
                    newparams = list(theta = c(0,0,0,1.0,0.5,0.3), beta = 2, sigma = 2))[[1]]

cat("range(y2):", range(dat$y2), "\n")

## build a devfun that takes c(theta, beta) jointly (nAGQ=1 style: only u
## is profiled internally by PIRLS, theta and beta are both externally
## supplied) -- this isolates a single PIRLS evaluation, decoupled from
## any particular optimizer's trajectory through parameter space
gm <- glFormula(form2, data = dat, family = Gamma(link = "log"))
devfun0 <- do.call(mkGlmerDevfun, gm)
devfun1 <- updateGlmerDevfun(devfun0, gm$reTrms, nAGQ = 1L)

theta_true <- c(0, 0, 0, 1.0, 0.5, 0.3)
beta_true <- 2

set.seed(42)
B <- 200
status <- character(B)
deviance <- rep(NA_real_, B)

for (i in seq_len(B)) {
  ## perturb by +/- 10%; elements that are exactly 0 in the true (singular)
  ## structure stay exactly 0 (0 +/- 10% of 0 == 0), which is correct: that
  ## reflects the deliberately-singular truth, not a bug in the sampling
  theta_i <- theta_true * (1 + runif(length(theta_true), -0.1, 0.1))
  beta_i  <- beta_true  * (1 + runif(1, -0.1, 0.1))

  val <- tryCatch(devfun1(c(theta_i, beta_i)), error = function(e) e)

  if (inherits(val, "error")) {
    status[i] <- "error"
  } else if (!is.finite(val)) {
    status[i] <- "bad fit (non-finite deviance)"
  } else {
    status[i] <- "good fit"
    deviance[i] <- val
  }
}

cat("\n=== results over", B, "perturbed (theta,beta) evaluations, +/-10% ===\n")
print(table(status))
cat("\ndeviance summary for 'good fit' evaluations:\n")
print(summary(deviance[status == "good fit"]))

outfile <- sprintf(
  "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/reprex_pirls_fragility_results_%s.rds",
  if (USE_OLD) "old" else "current")
saveRDS(list(status = status, deviance = deviance, theta_true = theta_true, beta_true = beta_true,
             version = as.character(packageVersion("lme4")), use_old = USE_OLD),
        outfile)
cat("\nsaved to", outfile, "\n")
