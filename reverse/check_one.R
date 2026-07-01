#!/usr/bin/env Rscript
##
## check_one.R -- run R CMD check on one reverse-dependency tarball
##
## Called by slurm_job.sh with SLURM_ARRAY_TASK_ID as the sole argument.
## The 1-based index selects a line from /opt/revdep/pkgs_to_check.txt.
##
## Output directory is named rdepends_PKGNAME.Rcheck so that the results
## are compatible with tools::check_packages_in_dir_changes() and the
## existing checkChanges.R script.
##
## Environment variables:
##   REVDEP_RESULTS  where to write .Rcheck dirs (default /results,
##                   bind-mounted by slurm_job.sh)
##   REVDEP_LME4     which lme4 version to use: "old" or "new" (default "new")
##                   selects /opt/revdep/Library_{old,new}/ as the first
##                   entry on .libPaths() so the right lme4 is found first

args <- commandArgs(trailingOnly = TRUE)
stopifnot("expected one argument: task index" = length(args) >= 1L)
idx <- as.integer(args[[1L]])
stopifnot(!is.na(idx), idx >= 1L)

## ---- Select lme4 version --------------------------------------------------
lme4_ver <- Sys.getenv("REVDEP_LME4", unset = "new")
stopifnot("REVDEP_LME4 must be 'old' or 'new'" = lme4_ver %in% c("old", "new"))
lme4_lib <- file.path("/opt/revdep", paste0("Library_", lme4_ver))
.libPaths(c(lme4_lib, .libPaths()))
cat(sprintf("lme4 version   : %s  (library: %s)\n", lme4_ver, lme4_lib))
cat(sprintf("lme4 installed : %s\n",
            as.character(packageVersion("lme4"))))

## ---- Identify the tarball -------------------------------------------------
PKG_LIST_FILE <- "/opt/revdep/pkgs_to_check.txt"
pkg_list      <- readLines(PKG_LIST_FILE)
stopifnot("task index out of range" = idx <= length(pkg_list))

tarball      <- pkg_list[[idx]]
pkg_name     <- sub("_.*$", "", basename(tarball))
results_dir  <- Sys.getenv("REVDEP_RESULTS", unset = "/results")
out_parent   <- results_dir
## rdepends_ prefix matches tools::check_packages_in_dir convention
out_name     <- sprintf("rdepends_%s.Rcheck", pkg_name)

dir.create(out_parent, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[task %d / %d]  package : %s\n", idx, length(pkg_list), pkg_name))
cat(sprintf("               tarball : %s\n", tarball))
cat(sprintf("               output  : %s\n", file.path(out_parent, out_name)))

## ---- Run R CMD check ------------------------------------------------------
## Set by slurm_job.sh and inherited by R CMD check:
##   _R_CHECK_FORCE_SUGGESTS_=false        : missing Suggests don't fail the check
##   _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE  : skip the network-dependent part of
##     the CRAN incoming feasibility check (compute nodes have no internet;
##     without this, every check hangs/fails at "checking CRAN incoming
##     feasibility")

## R CMD check places PKGNAME.Rcheck under the directory given by -o;
## rename to rdepends_PKGNAME.Rcheck afterwards for checkChanges.R compat.
ret <- system2("R", c("CMD", "check", "--as-cran", "--no-manual",
                       "-o", out_parent, tarball))

std_path  <- file.path(out_parent, paste0(pkg_name, ".Rcheck"))
rdep_path <- file.path(out_parent, out_name)
if (file.exists(rdep_path))
    unlink(rdep_path, recursive = TRUE)  # remove stale results from a previous run
if (file.exists(std_path))
    file.rename(std_path, rdep_path)

## Always exit 0 so SLURM marks the task COMPLETED regardless of check
## outcome; R CMD check results (ERROR/WARNING/NOTE) live in the .Rcheck
## directory and are evaluated by checkChanges.R, not by SLURM exit status.
cat(sprintf("R CMD check exit code for %s: %d\n", pkg_name, ret))
quit(status = 0L)
