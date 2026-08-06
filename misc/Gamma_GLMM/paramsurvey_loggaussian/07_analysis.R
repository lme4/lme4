## Combine the five per-method result files and summarize: status/singular
## counts, timing, and parameter recovery (beta, RE sd, sigma, -2*logLik)
## against the "pretty" true values used to simulate the data. Note: PQL's
## negll is always NA (not a real ML method -- see toolkit.R).

wd <- "/home/bolker/R/pkgs/lme4/misc/Gamma_GLMM/paramsurvey_loggaussian"
sim <- readRDS(file.path(wd, "rail_simdata.rds"))

methods <- c("glmmTMB", "jointphi", "glmer", "mgcv", "pql")
res <- setNames(lapply(methods, function(m) {
  f <- file.path(wd, paste0("rail_results_", m, ".rds"))
  if (!file.exists(f)) return(NULL)
  readRDS(f)
}), methods)
res <- res[!vapply(res, is.null, logical(1))]

cat("=== true (pretty) parameters ===\n")
cat(sprintf("beta = %.4f  sd = %.4f  sigma = %.4f\n",
            sim$pretty$beta, sim$pretty$sd, sim$pretty$sigma))

cat("\n=== status counts ===\n")
for (m in names(res)) { cat(m, ": "); print(table(res[[m]]$status)) }

cat("\n=== singular counts ===\n")
for (m in names(res)) { cat(m, ": "); print(table(res[[m]]$singular)) }

cat("\n=== median elapsed time (s) ===\n")
for (m in names(res)) cat(sprintf("%-10s %.4f\n", m, median(res[[m]]$time_sec, na.rm = TRUE)))

cat("\n=== parameter recovery (median over replicates; true value in parens) ===\n")
for (m in names(res)) {
  d <- res[[m]]
  cat(sprintf("\n%s:\n", m))
  cat(sprintf("  beta  = %.4f (%.2f)\n", median(d$beta, na.rm = TRUE), sim$pretty$beta))
  cat(sprintf("  sd    = %.4f (%.2f)\n", median(d$sd, na.rm = TRUE), sim$pretty$sd))
  cat(sprintf("  sigma = %.4f (%.2f)\n", median(d$sigma, na.rm = TRUE), sim$pretty$sigma))
  cat(sprintf("  negll = %.4f\n", median(d$negll, na.rm = TRUE)))
}

cat("\n=== paired per-replicate sigma, all methods (wide) ===\n")
wide <- data.frame(i = seq_len(sim$B))
for (m in names(res)) wide[[m]] <- res[[m]]$sigma[match(wide$i, res[[m]]$i)]
print(wide)

saveRDS(res, file.path(wd, "rail_results_combined.rds"))
cat("\nsaved combined results to", file.path(wd, "rail_results_combined.rds"), "\n")
