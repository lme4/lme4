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
## Environment variable REVDEP_RESULTS controls where results land
## (default /results, which is bind-mounted by slurm_job.sh).

args <- commandArgs(trailingOnly = TRUE)
stopifnot("expected one argument: task index" = length(args) >= 1L)
idx <- as.integer(args[[1L]])
stopifnot(!is.na(idx), idx >= 1L)

PKG_LIST_FILE <- "/opt/revdep/pkgs_to_check.txt"
pkg_list      <- readLines(PKG_LIST_FILE)
stopifnot("task index out of range" = idx <= length(pkg_list))

tarball      <- pkg_list[[idx]]
pkg_name     <- sub("_.*$", "", basename(tarball))
results_dir  <- Sys.getenv("REVDEP_RESULTS", unset = "/results")
## use rdepends_ prefix to match tools::check_packages_in_dir convention
out_parent   <- results_dir
out_name     <- sprintf("rdepends_%s.Rcheck", pkg_name)

dir.create(out_parent, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[task %d / %d]  package : %s\n", idx, length(pkg_list), pkg_name))
cat(sprintf("               tarball : %s\n", tarball))
cat(sprintf("               output  : %s\n", file.path(out_parent, out_name)))

## Match the check environment used in checkReverse.R
old_sug <- Sys.getenv("_R_CHECK_FORCE_SUGGESTS_", unset = NA_character_)
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = "false")
on.exit({
    if (is.na(old_sug)) Sys.unsetenv("_R_CHECK_FORCE_SUGGESTS_")
    else Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = old_sug)
}, add = TRUE)

## R CMD check places PKGNAME.Rcheck under the directory given by -o;
## rename to rdepends_PKGNAME.Rcheck afterwards for checkChanges.R compat.
ret <- system2("R", c("CMD", "check", "--as-cran", "--no-manual",
                       "-o", out_parent, tarball))

std_path  <- file.path(out_parent, paste0(pkg_name, ".Rcheck"))
rdep_path <- file.path(out_parent, out_name)
if (file.exists(std_path))
    file.rename(std_path, rdep_path)

cat(sprintf("exit code for %s: %d\n", pkg_name, ret))
quit(status = ret)
