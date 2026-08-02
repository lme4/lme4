## Fit lme4 2.0-6 (unmodified, isolated install) to each simulated dataset
## via plain glmer(). Must load the old lme4 BEFORE sourcing toolkit.R,
## since toolkit.R's own library(lme4) call is a no-op once a package of
## that name is already attached (R doesn't reload from a different
## lib.loc for an already-attached package).

tmplib <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/lme4_206_lib"
library(lme4, lib.loc = tmplib)
stopifnot(packageVersion("lme4") == "2.0-6")

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_complex"
MC_CORES <- if (length(args) >= 2) as.integer(args[[2]]) else 1

wd <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"
source(file.path(wd, "toolkit.R"))
library(parallel)
cat("Using lme4", as.character(packageVersion("lme4")), "(unmodified, from isolated lib) -- mc.cores =", MC_CORES, "\n")

sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
cat("Fitting lme4 2.0-6 to", sim$B, "replicates of", sim$name, "\n")

applyfun <- if (MC_CORES > 1) function(...) mclapply(..., mc.cores = MC_CORES) else lapply
results <- applyfun(seq_len(sim$B), function(b) fit_glmer_one(b, sim$formula, sim$sim_data[[b]], sim$family))
df <- resultsToDF(results, names(sim$pretty$beta))

cat("\n=== status ===\n"); print(table(df$status))
cat("=== singular ===\n"); print(table(df$singular))
print(df[, c("i", "status", "singular", "time_sec", "sd1", "sd2", "corr", "phi", "negll")])

saveRDS(df, file.path(wd, paste0(example, "_results_lme4old.rds")))
cat("\nsaved to", file.path(wd, paste0(example, "_results_lme4old.rds")), "\n")
