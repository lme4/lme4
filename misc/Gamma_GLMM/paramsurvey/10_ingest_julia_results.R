## Read the Julia-side (MixedModels.jl, pa/dispersion-again) fit results
## CSV and reshape into the same resultsToDF() shape used by the R-side
## methods (02-06), so 07_analysis.R can combine all methods uninterrupted.
##
## Expects julia/results_<example>.csv with columns: i, status, singular,
## msg, time_sec, sd1, sd2, corr, phi, negll, plus one column per fixed
## effect IN THE SAME ORDER as sim$pretty$beta (Julia's coefnames aren't
## forced to match R's factor-contrast column names character-for-character,
## so we match fixed effects by position, not by name).

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_simple"

wd <- here::here("misc/Gamma_GLMM/paramsurvey")
sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
beta_names <- names(sim$pretty$beta)

infile <- file.path(wd, "julia", paste0("results_", example, ".csv"))
raw <- read.csv(infile, stringsAsFactors = FALSE)

stopifnot(nrow(raw) == sim$B)
p <- length(beta_names)
beta_cols <- paste0("beta", seq_len(p))
stopifnot(all(beta_cols %in% names(raw)))

df <- data.frame(i = raw$i, status = raw$status, singular = as.logical(raw$singular),
                  msg = raw$msg, time_sec = raw$time_sec, sd1 = raw$sd1, sd2 = raw$sd2,
                  corr = raw$corr, phi = raw$phi, negll = raw$negll,
                  raw[, beta_cols, drop = FALSE], check.names = FALSE)
names(df)[names(df) %in% beta_cols] <- beta_names

cat("\n=== status ===\n"); print(table(df$status))
cat("=== singular ===\n"); print(table(df$singular))
print(df[, c("i", "status", "singular", "time_sec", "sd1", "sd2", "corr", "phi", "negll")])

saveRDS(df, file.path(wd, paste0(example, "_results_juliaMixedModels.rds")))
cat("\nsaved to", file.path(wd, paste0(example, "_results_juliaMixedModels.rds")), "\n")
