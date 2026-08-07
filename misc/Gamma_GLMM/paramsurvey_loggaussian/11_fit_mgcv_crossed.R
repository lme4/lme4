wd <- here::here("misc", "Gamma_GLMM", "paramsurvey_loggaussian")
source(file.path(wd, "toolkit.R"))

family_name <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else "gaussian"
sim <- readRDS(file.path(wd, sprintf("crossed_simdata_%s.rds", family_name)))

for (case in names(sim$cases)) {
  cat("=== case:", case, "===\n")
  sim_data <- sim$cases[[case]]$sim_data
  fitfun <- if (grepl("^oneway", case)) fit_mgcv_oneway_one else fit_mgcv_crossed_one
  results <- lapply(seq_len(sim$B), function(b) fitfun(b, sim_data[[b]], family = sim$family))
  df <- resultsToDFCrossed(results)

  cat("\n=== status ===\n"); print(table(df$status))
  cat("=== singular ===\n"); print(table(df$singular))
  print(df[, c("i", "status", "singular", "time_sec", "beta", "sd1", "sd2", "sigma", "negll", "q")])

  outfile <- file.path(wd, sprintf("crossed_results_mgcv_%s_%s.rds", case, family_name))
  saveRDS(df, outfile)
  cat("\nsaved to", outfile, "\n\n")
}
