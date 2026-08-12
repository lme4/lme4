## Export one example's simulated replicate datasets to CSV so the Julia
## side (MixedModels.jl, pa/dispersion-again) can fit them without going
## through RCall/RDS. Only the columns the model formula actually needs
## are written (y, and whatever fixed/random-effect grouping columns
## appear in the formula).

args <- commandArgs(trailingOnly = TRUE)
example <- if (length(args) >= 1) args[[1]] else "epil2_simple"

wd <- here::here("misc/Gamma_GLMM/paramsurvey")
sim <- readRDS(file.path(wd, paste0(example, "_simdata.rds")))

outdir <- file.path(wd, "julia", "data", example)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vars <- all.vars(sim$formula)
cat("Exporting", sim$B, "replicates of", sim$name, "( columns:", vars, ") to", outdir, "\n")

for (b in seq_len(sim$B)) {
  write.csv(sim$sim_data[[b]][, vars], file.path(outdir, sprintf("rep_%03d.csv", b)),
            row.names = FALSE)
}

cat("done.\n")
