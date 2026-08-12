## Read one SAS-side (PROC GLIMMIX) results CSV and reshape into the
## same resultsToDF() shape the other methods use, so 07_analysis.R can
## combine it uninterrupted. Mirrors 10_ingest_julia_results.R exactly,
## including the positional (not name-matched) beta handling -- SAS's
## own effect-name strings (e.g. "trt" for R's "trtprogabide") don't
## string-match R's factor-contrast column names either.
##
## Expects sas/results_<example>_<method>.csv (method = rspl | laplace)
## with columns: i, status, singular, msg, time_sec, sd1, sd2, corr, phi,
## negll, beta1..betaN. Saves as <example>_results_sas<Rspl|Laplace>.rds.

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_simple"
method <- if (length(args) >= 2) args[[2]] else "rspl"
stopifnot(method %in% c("rspl", "laplace"))
out_method <- if (method == "rspl") "sasRSPL" else "sasLaplace"

wd <- here::here("misc/Gamma_GLMM/paramsurvey")
sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))
beta_names <- names(sim$pretty$beta)

infile <- file.path(wd, "sas", paste0("results_", example, "_", method, ".csv"))
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

outfile <- file.path(wd, paste0(example, "_results_", out_method, ".rds"))
saveRDS(df, outfile)
cat("\nsaved to", outfile, "\n")
