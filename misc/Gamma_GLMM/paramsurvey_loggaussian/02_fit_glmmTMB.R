wd <- "/home/bolker/R/pkgs/lme4/misc/Gamma_GLMM/paramsurvey_loggaussian"
source(file.path(wd, "toolkit.R"))

sim <- readRDS(file.path(wd, "rail_simdata.rds"))
cat("Fitting glmmTMB to", sim$B, "replicates of", sim$name, "\n")

results <- lapply(seq_len(sim$B), function(b) fit_glmmTMB_one(b, sim$sim_data[[b]]))
df <- resultsToDF(results)

cat("\n=== status ===\n"); print(table(df$status))
cat("=== singular ===\n"); print(table(df$singular))
print(df[, c("i", "status", "singular", "time_sec", "beta", "sd", "sigma", "negll")])

saveRDS(df, file.path(wd, "rail_results_glmmTMB.rds"))
cat("\nsaved to", file.path(wd, "rail_results_glmmTMB.rds"), "\n")
