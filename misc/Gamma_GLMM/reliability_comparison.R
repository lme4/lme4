## Compare the proportion of glmer() fits that error / warn / converge
## cleanly across many *similar* (not identical) replicate datasets, for
## unmodified lme4 2.0-6 vs. the current (dev, Gamma_GLMM branch) lme4.
## Each replicate resimulates the whole dataset (design + response) fresh
## from the same generative process used in tests/testthat/test-isSingular.R
## ("checking singular fit for merMod"): 20x20 balanced design, two crossed
## random-slope terms, one deliberately singular (group1: theta=0,0,0),
## one not (group2: theta=1.0,0.5,0.3), Gamma(link=log) response with
## shape=0.25 (sigma=2).
##
## Uses a shared seed sequence so both library versions see the *same* B
## datasets (a fair paired comparison), run as two separate processes (each
## R session can only have one lme4 loaded) and combined afterward.

USE_OLD <- FALSE   # TRUE: lme4 2.0-6 (unmodified); FALSE: current dev install

if (USE_OLD) {
  tmplib <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/lme4_206_lib"
  library(lme4, lib.loc = tmplib)
  stopifnot(packageVersion("lme4") == "2.0-6")
} else {
  library(lme4)
}
cat("Using lme4", as.character(packageVersion("lme4")),
    if (USE_OLD) "(unmodified, from isolated lib)" else "(current dev install)", "\n")

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

B <- 50
master_seed <- 9000
status <- character(B)
msgs <- vector("list", B)

for (b in seq_len(B)) {
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
  msgs[[b]] <- warn_msgs
  if (inherits(fit, "error")) {
    status[b] <- "error"
  } else if (length(warn_msgs) > 0) {
    status[b] <- "converged with warnings"
  } else {
    status[b] <- "clean convergence"
  }
  cat(sprintf("rep %d: %s\n", b, status[b]))
}

cat("\n=== summary over", B, "replicate datasets ===\n")
print(table(status))

outfile <- sprintf(
  "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/report4bb/reliability_comparison_%s.rds",
  if (USE_OLD) "old" else "current")
saveRDS(list(status = status, msgs = msgs, version = as.character(packageVersion("lme4")), use_old = USE_OLD, B = B),
        outfile)
cat("\nsaved to", outfile, "\n")
