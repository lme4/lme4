## Usage: Rscript run_dispsim.R [nsim] [ngrp] [nper] [shape] [sd_re] [workers]
## Run from misc/Gamma_GLMM/dispsim/. Writes dispsim_<setup>.rds.

## one BLAS/OpenMP thread per worker (also inherited by forked workers)
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

source("dispsim_funs.R")
library(furrr)

args <- commandArgs(trailingOnly = TRUE)
argn <- function(k, default) if (length(args) >= k) as.numeric(args[k]) else default
nsim    <- argn(1, 100)
ngrp    <- argn(2, 20)
nper    <- argn(3, 20)
shape   <- argn(4, 0.5)
sd_re   <- argn(5, 0.5)
workers <- argn(6, 12)

plan(multicore, workers = workers)
setup <- sprintf("g%d_n%d_shape%g_sd%g", ngrp, nper, shape, sd_re)
t0 <- Sys.time()
res <- future_map_dfr(seq_len(nsim), run_one,
                      ngrp = ngrp, nper = nper, shape = shape, sd_re = sd_re,
                      .options = furrr_options(seed = NULL))
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
attr(res, "truth") <- c(beta0 = 1, beta1 = 0.5, sd_re = sd_re, phi = 1 / shape)
attr(res, "setup") <- list(nsim = nsim, ngrp = ngrp, nper = nper, shape = shape,
                           sd_re = sd_re, workers = workers, wall = wall,
                           lme4 = as.character(packageVersion("lme4")),
                           glmmTMB = as.character(packageVersion("glmmTMB")))
saveRDS(res, sprintf("dispsim_%s.rds", setup))
cat(sprintf("%s: %d sims, %.1f s wall on %d workers\n", setup, nsim, wall, workers))
