## Combine the four per-method result files for one example and summarize:
## status/singular counts, timing, and parameter recovery (RE sd/corr,
## dispersion, fixed effects) against the "pretty" true values used to
## simulate the data.

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_complex"

wd <- "/tmp/claude-1000/-home-bolker-Documents-R-pkgs-lme4/6bffb877-f2f1-42e7-974f-8b4550ca1ead/scratchpad/param_survey"
sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))

methods <- c("glmmTMB", "jointphi", "pirlsdigamma", "pirlsmoment", "lme4current", "lme4old")
res <- setNames(lapply(methods, function(m) {
  f <- file.path(wd, paste0(example, "_results_", m, ".rds"))
  if (!file.exists(f)) return(NULL)
  readRDS(f)
}), methods)
res <- res[!vapply(res, is.null, logical(1))]

cat("=== true (pretty) parameters ===\n")
cat("sd1 =", sim$pretty$sd[1], " sd2 =", sim$pretty$sd[2],
    " corr =", sim$pretty$corr, " phi =", sim$pretty$phi, "\n")
print(sim$pretty$beta)

cat("\n=== status counts ===\n")
for (m in names(res)) { cat(m, ": "); print(table(res[[m]]$status)) }

cat("\n=== singular counts ===\n")
for (m in names(res)) { cat(m, ": "); print(table(res[[m]]$singular)) }

cat("\n=== median elapsed time (s) ===\n")
for (m in names(res)) cat(sprintf("%-12s %.3f\n", m, median(res[[m]]$time_sec, na.rm = TRUE)))

cat("\n=== parameter recovery (median over replicates; true value in parens) ===\n")
for (m in names(res)) {
  d <- res[[m]]
  cat(sprintf("\n%s:\n", m))
  cat(sprintf("  sd1  = %.4f (%.2f)\n", median(d$sd1, na.rm = TRUE), sim$pretty$sd[1]))
  cat(sprintf("  sd2  = %.4f (%.2f)\n", median(d$sd2, na.rm = TRUE), sim$pretty$sd[2]))
  cat(sprintf("  corr = %.4f (%.2f)\n", median(d$corr, na.rm = TRUE), sim$pretty$corr))
  cat(sprintf("  phi  = %.4f (%.2f)\n", median(d$phi, na.rm = TRUE), sim$pretty$phi))
}

saveRDS(res, file.path(wd, paste0(example, "_results_combined.rds")))
cat("\nsaved combined results to", file.path(wd, paste0(example, "_results_combined.rds")), "\n")
