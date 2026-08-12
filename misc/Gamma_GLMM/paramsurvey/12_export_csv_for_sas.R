## Export one example's simulated replicate datasets to a single stacked
## CSV (all B replicates, keyed by a `rep` column) so the SAS side
## (PROC GLIMMIX) can fit all of them in one BY-group pass rather than
## one file/invocation per replicate -- SAS's natural idiom for "fit the
## same model to many groups" is `BY rep;`, not a loop over files the
## way the Julia side's DaemonMode client does it. Only the columns the
## model formula actually needs are written.

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_simple"

wd <- here::here("misc/Gamma_GLMM/paramsurvey")
sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))

outdir <- file.path(wd, "sas", "data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vars <- all.vars(sim$formula)
cat("Exporting", sim$B, "replicates of", sim$name, "( columns: rep,", vars, ") to", outdir, "\n")

stacked <- do.call(rbind, lapply(seq_len(sim$B), function(b) {
  d <- sim$sim_data[[b]][, vars]
  cbind(rep = b, d)
}))

outfile <- file.path(outdir, paste0(example, ".csv"))
write.csv(stacked, outfile, row.names = FALSE)
cat("done:", nrow(stacked), "rows ->", outfile, "\n")
