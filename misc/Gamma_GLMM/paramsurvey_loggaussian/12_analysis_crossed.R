## Combine the crossed-RE per-method, per-case result files and summarize:
## status/singular counts, timing, parameter recovery (beta, sd1, sd2,
## sigma, -2*logLik) against the true simulation parameters, and a
## cross-method sanity check that q (total RE levels) agrees for glmmTMB/
## glmer/mgcv on the same replicate (it's a property of the data, so any
## mismatch would indicate a method silently dropping levels).

wd <- here::here("misc", "Gamma_GLMM", "paramsurvey_loggaussian")
family_name <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else "gaussian"
sim <- readRDS(file.path(wd, sprintf("crossed_simdata_%s.rds", family_name)))

methods <- c("glmmTMB", "glmer", "mgcv")
cases <- names(sim$cases)

res <- setNames(lapply(cases, function(cs) {
  setNames(lapply(methods, function(m) {
    f <- file.path(wd, sprintf("crossed_results_%s_%s_%s.rds", m, cs, family_name))
    if (!file.exists(f)) return(NULL)
    readRDS(f)
  }), methods)
}), cases)

cat("=== true (pretty) parameters ===\n")
cat(sprintf("beta = %.4f  sd1 = %.4f  sd2 = %.4f  sigma = %.4f\n",
            sim$pretty$beta, sim$pretty$sd1, sim$pretty$sd2, sim$pretty$sigma))

for (cs in cases) {
  cat(sprintf("\n\n########## case: %s ##########\n", cs))
  r <- res[[cs]][!vapply(res[[cs]], is.null, logical(1))]

  cat("\n=== status counts ===\n")
  for (m in names(r)) { cat(m, ": "); print(table(r[[m]]$status)) }

  cat("\n=== singular counts ===\n")
  for (m in names(r)) { cat(m, ": "); print(table(r[[m]]$singular)) }

  cat("\n=== median elapsed time (s) ===\n")
  for (m in names(r)) cat(sprintf("%-10s %.4f\n", m, median(r[[m]]$time_sec, na.rm = TRUE)))

  cat("\n=== parameter recovery (median over replicates; true value in parens) ===\n")
  for (m in names(r)) {
    d <- r[[m]]
    cat(sprintf("\n%s:\n", m))
    cat(sprintf("  beta  = %.4f (%.2f)\n", median(d$beta, na.rm = TRUE), sim$pretty$beta))
    cat(sprintf("  sd1   = %.4f (%.2f)\n", median(d$sd1, na.rm = TRUE), sim$pretty$sd1))
    cat(sprintf("  sd2   = %.4f (%.2f)\n", median(d$sd2, na.rm = TRUE), sim$pretty$sd2))
    cat(sprintf("  sigma = %.4f (%.2f)\n", median(d$sigma, na.rm = TRUE), sim$pretty$sigma))
    if (!all(is.na(d$corr)))
      cat(sprintf("  corr  = %.4f (%.2f)\n", median(d$corr, na.rm = TRUE), sim$pretty$corr))
    cat(sprintf("  negll = %.4f\n", median(d$negll, na.rm = TRUE)))
  }

  cat("\n=== q (total RE levels) agreement across methods (should be identical per replicate) ===\n")
  wide_q <- data.frame(i = seq_len(sim$B))
  for (m in names(r)) wide_q[[m]] <- r[[m]]$q[match(wide_q$i, r[[m]]$i)]
  print(wide_q)

  cat("\n=== paired per-replicate sigma, all methods (wide) ===\n")
  wide <- data.frame(i = seq_len(sim$B))
  for (m in names(r)) wide[[m]] <- r[[m]]$sigma[match(wide$i, r[[m]]$i)]
  print(wide)
}

outfile <- file.path(wd, sprintf("crossed_results_combined_%s.rds", family_name))
saveRDS(res, outfile)
cat("\nsaved combined results to", outfile, "\n")
