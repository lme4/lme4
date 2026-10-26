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
##   REVDEP_KEEP_ALL "true" keeps the whole .Rcheck directory instead of
##                   pruning it to the logs (default "false"; see below)

args <- commandArgs(trailingOnly = TRUE)
stopifnot("expected one argument: task index" = length(args) >= 1L)
idx <- as.integer(args[[1L]])
stopifnot(!is.na(idx), idx >= 1L)

## ---- Select lme4 version --------------------------------------------------
lme4_ver <- Sys.getenv("REVDEP_LME4", unset = "new")
stopifnot("REVDEP_LME4 must be 'old' or 'new'" = lme4_ver %in% c("old", "new"))
lme4_lib <- file.path("/opt/revdep", paste0("Library_", lme4_ver))
.libPaths(c(lme4_lib, .libPaths()))
## R CMD check runs in child R processes (system2() below), which do not
## inherit .libPaths(); without R_LIBS they find the r2u-installed lme4 in
## /usr/lib/R/site-library (pulled in as an apt dependency) instead, so
## 'old' and 'new' would both silently check against that version.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
cat(sprintf("lme4 version   : %s  (library: %s)\n", lme4_ver, lme4_lib))
cat(sprintf("lme4 installed : %s\n",
            as.character(packageVersion("lme4"))))

## ---- Per-package user cache -----------------------------------------------
## Some packages (e.g. RBioFormats) download files into
## tools::R_user_dir(pkg, "cache") on first load. Under singularity --no-home
## that is not a usable location and the image itself is read-only, so link
## the copies pre-downloaded at build time (/opt/revdep/Rcache) into a
## writable per-job directory; R CMD check inherits R_USER_CACHE_DIR.
cache_seed <- "/opt/revdep/Rcache"
cache_dir  <- file.path(tempdir(), "Rcache")
dir.create(cache_dir, showWarnings = FALSE)
for (f in list.files(cache_seed, recursive = TRUE)) {
    dir.create(dirname(file.path(cache_dir, f)), recursive = TRUE,
               showWarnings = FALSE)
    file.symlink(file.path(cache_seed, f), file.path(cache_dir, f))
}
Sys.setenv(R_USER_CACHE_DIR = cache_dir)

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

## ---- Prune the .Rcheck directory ------------------------------------------
## A full .Rcheck directory (unpacked sources, installed package, vignette
## and test outputs) can hold thousands of files; a full run of ~650
## packages x old/new then exhausts the file-count quota on /scratch, and
## once it does, SLURM can no longer create the tasks' log files and every
## remaining task fails (exit code 0:53). checkChanges.R,
## checkChanges_detailed.R and grep_problems.sh only read 00check.log and
## 00install.out; example output and failed test output are kept too, for
## diagnosing failures. Set REVDEP_KEEP_ALL=true to keep everything.
keep_all <- tolower(Sys.getenv("REVDEP_KEEP_ALL", "false")) %in% c("true", "1", "yes")
if (!keep_all && dir.exists(rdep_path)) {
    all_files <- list.files(rdep_path, recursive = TRUE, all.files = TRUE)
    keep_re   <- "^00check\\.log$|^00install\\.out$|-Ex\\.Rout$|\\.Rout\\.fail$"
    keep      <- all_files[grepl(keep_re, all_files) &
                           !startsWith(all_files, "00_pkg_src/")]
    save_dir  <- tempfile("keep")
    for (f in keep) {
        dir.create(dirname(file.path(save_dir, f)), recursive = TRUE,
                   showWarnings = FALSE)
        file.copy(file.path(rdep_path, f), file.path(save_dir, f))
    }
    unlink(rdep_path, recursive = TRUE)
    dir.create(rdep_path)
    for (f in keep) {
        dir.create(dirname(file.path(rdep_path, f)), recursive = TRUE,
                   showWarnings = FALSE)
        file.copy(file.path(save_dir, f), file.path(rdep_path, f))
    }
    unlink(save_dir, recursive = TRUE)
    cat(sprintf("pruned %s: kept %d of %d files\n",
                out_name, length(keep), length(all_files)))
}

## Always exit 0 so SLURM marks the task COMPLETED regardless of check
## outcome; R CMD check results (ERROR/WARNING/NOTE) live in the .Rcheck
## directory and are evaluated by checkChanges.R, not by SLURM exit status.
cat(sprintf("R CMD check exit code for %s: %d\n", pkg_name, ret))
quit(status = 0L)
