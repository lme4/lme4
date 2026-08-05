## Fit the R-level nested-fixed-point PIRLS devfun (toolkit.R's
## make_pirls_phi_devfun) to each simulated dataset, for a given phiType:
##   digamma -- true Gamma conditional MLE for phi each outer iteration
##              (README #6's "ALSO FALSIFIED" idea, resurrected as arm 5
##              to check whether that single-RE-scenario null result
##              generalizes to the richer models in this survey)
##   moment  -- the same crude dev/n plug-in the C++ "current" fix uses,
##              but implemented in R instead of C++ -- comparing this
##              against the lme4current arm (identical algorithm, same
##              starting values/bounds where applicable) isolates the pure
##              R-vs-compiled-C++ computational overhead, since any
##              parameter-recovery/reliability differences between them
##              should be ~0 (same math) and only speed should differ.

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_complex"
phiType <- if (length(args) >= 2) args[[2]] else "digamma"
stopifnot(phiType %in% c("digamma", "moment"))
MC_CORES <- if (length(args) >= 3) as.integer(args[[3]]) else 1
## maxPhiIter: cap on the outer nested-fixed-point loop (default 30, i.e.
## "run to convergence"). Set to 1 to test the single-update-per-devfun-call
## variant (TODO.md #2: is it worth skipping/cheapening the nested loop?)
## -- phi is still incorporated correctly into the PIRLS deviance
## expression, it just isn't iterated to convergence against u/beta within
## one devfun evaluation.
maxPhiIter <- if (length(args) >= 4) as.integer(args[[4]]) else 30L

wd <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"
source(file.path(wd, "toolkit.R"))
library(parallel)

sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
spec <- list(formula = sim$formula, family = sim$family)
cat("Fitting PIRLS/", phiType, "-phi (maxPhiIter=", maxPhiIter, ") to ", sim$B,
    " replicates of ", sim$name, " -- mc.cores = ", MC_CORES, "\n", sep = "")

applyfun <- if (MC_CORES > 1) function(...) mclapply(..., mc.cores = MC_CORES) else lapply
results <- applyfun(seq_len(sim$B), function(b) fit_pirls_phi_one(b, spec, sim$sim_data[[b]], phiType = phiType, maxPhiIter = maxPhiIter))
df <- resultsToDF(results, names(sim$pretty$beta))

cat("\n=== status ===\n"); print(table(df$status))
cat("=== singular ===\n"); print(table(df$singular))
print(df[, c("i", "status", "singular", "time_sec", "sd1", "sd2", "corr", "phi", "negll")])

label <- if (phiType == "digamma") "pirlsdigamma" else "pirlsmoment"
if (maxPhiIter != 30L) label <- paste0(label, "_maxPhiIter", maxPhiIter)
saveRDS(df, file.path(wd, paste0(example, "_results_", label, ".rds")))
cat("\nsaved to", file.path(wd, paste0(example, "_results_", label, ".rds")), "\n")
