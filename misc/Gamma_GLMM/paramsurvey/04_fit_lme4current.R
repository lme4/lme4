## Fit lme4 current (dev, moment-phi C++ fix, main env install) to each
## simulated dataset via plain glmer().

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_complex"
MC_CORES <- if (length(args) >= 2) as.integer(args[[2]]) else 1
## maxPhiIter: override glmerControl()'s default (100) for the real C++
## nested-fixed-point phi loop. Leave unset (NA) to use glmer()'s own
## default control, preserving the existing baseline output file.
maxPhiIter <- if (length(args) >= 3) as.integer(args[[3]]) else NA_integer_

wd <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"
source(file.path(wd, "toolkit.R"))
library(parallel)
cat("Using lme4", as.character(packageVersion("lme4")), "(current dev install) -- mc.cores =", MC_CORES,
    " maxPhiIter =", if (is.na(maxPhiIter)) "default" else maxPhiIter, "\n")

control <- if (is.na(maxPhiIter)) NULL else glmerControl(maxPhiIter = maxPhiIter)

sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
cat("Fitting lme4 current to", sim$B, "replicates of", sim$name, "\n")

applyfun <- if (MC_CORES > 1) function(...) mclapply(..., mc.cores = MC_CORES) else lapply
results <- applyfun(seq_len(sim$B), function(b) fit_glmer_one(b, sim$formula, sim$sim_data[[b]], sim$family, control = control))
df <- resultsToDF(results, names(sim$pretty$beta))

cat("\n=== status ===\n"); print(table(df$status))
cat("=== singular ===\n"); print(table(df$singular))
print(df[, c("i", "status", "singular", "time_sec", "sd1", "sd2", "corr", "phi", "negll")])

label <- if (is.na(maxPhiIter)) "lme4current" else paste0("lme4current_maxPhiIter", maxPhiIter)
saveRDS(df, file.path(wd, paste0(example, "_results_", label, ".rds")))
cat("\nsaved to", file.path(wd, paste0(example, "_results_", label, ".rds")), "\n")
