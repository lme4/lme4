#!/usr/bin/env Rscript
##
## setup_revdeps.R -- internet-connected setup phase, run during 'docker build'
##
## Adapted from checkReverse.R (the install-only / --no-check phase).
##
## Expects exactly two lme4_*.tar.gz files in /opt/revdep/ (old and new
## versions to compare).  Determines which is older by version number.
##
## Steps:
##   1. Configure CRAN + Bioconductor repositories
##   2. Discover all direct reverse dependencies of lme4
##   3. Download their source tarballs into /opt/revdep/tarballs/
##   4. Install all transitive dependencies into the container R library
##      (lme4 itself is excluded here)
##   5. Install old lme4 into /opt/revdep/Library_old/
##   6. Install new lme4 into /opt/revdep/Library_new/

options(timeout = 120)   # default 60s can time out on slow connections

REVDEP_DIR    <- "/opt/revdep"
TARBALL_DIR   <- file.path(REVDEP_DIR, "tarballs")
PKG_LIST_FILE <- file.path(REVDEP_DIR, "pkgs_to_check.txt")
LIB_OLD       <- file.path(REVDEP_DIR, "Library_old")
LIB_NEW       <- file.path(REVDEP_DIR, "Library_new")
NCPUS         <- max(1L, parallel::detectCores())

for (d in c(TARBALL_DIR, LIB_OLD, LIB_NEW))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- Locate and sort the two lme4 tarballs --------------------------------
lme4_tgz <- Sys.glob(file.path(REVDEP_DIR, "lme4_*.tar.gz"))
stopifnot("expected exactly two lme4 tarballs (old and new)" =
              length(lme4_tgz) == 2L)
lme4_vers <- package_version(
    sub("^lme4_(.*)\\.tar\\.gz$", "\\1", basename(lme4_tgz)))
ord      <- order(lme4_vers)
lme4_old <- lme4_tgz[ord[1L]]
lme4_new <- lme4_tgz[ord[2L]]
cat(sprintf("old lme4: %s\nnew lme4: %s\n\n",
            basename(lme4_old), basename(lme4_new)))

## ---- 1. Repos -------------------------------------------------------------
with_bioc <- identical(Sys.getenv("WITH_BIOC", "true"), "true")
if (with_bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    repos <- BiocManager::repositories()
} else {
    repos <- c(CRAN = "https://cloud.r-project.org")
}
options(repos = repos,
        Ncpus = NCPUS,
        install.packages.compile.from.source = "always",
        useFancyQuotes = FALSE)
cat(sprintf("Bioconductor   : %s\n", if (with_bioc) "included" else "skipped"))
cat("Repositories:\n"); print(repos); cat("\n")

## ---- 2. Discover reverse dependencies of lme4 ----------------------------
cat("--- Discovering reverse dependencies of lme4 ---\n")
ap <- available.packages(repos = repos)
rdeps <- sort(tools::package_dependencies(
    "lme4", db = ap,
    which   = c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances"),
    reverse = TRUE)[[1L]])
cat(sprintf("Found %d reverse dependencies\n", length(rdeps)))

## ---- 3. Download source tarballs -----------------------------------------
## download.packages() silently skips packages not found in any repo and
## returns only successfully located packages, so we detect failures by
## comparing the returned names against rdeps.
cat("\n--- Downloading source tarballs ---\n")
dl <- download.packages(rdeps, destdir = TARBALL_DIR,
                        repos = repos, type = "source")
failed <- setdiff(rdeps, dl[, 1L])
if (length(failed))
    warning(sprintf("Could not download %d package(s): %s",
                    length(failed), paste(failed, collapse = ", ")))
cat(sprintf("Downloaded %d / %d tarballs to %s\n",
            nrow(dl), length(rdeps), TARBALL_DIR))
writeLines(dl[, 2L], PKG_LIST_FILE)

## ---- 4. Install all transitive dependencies (excluding lme4) -------------
## Resolve the full dependency closure of the downloaded rev-deps so that
## R CMD check can run offline on Compute Canada.
## r2u installs binary apt packages where available (fast); falls back to source.
##
## Recursive resolution uses only hard deps (Depends/Imports/LinkingTo) to
## avoid pulling in the enormous transitive closure of all Suggests.
## Suggests of the rev-dep packages themselves are included one level deep
## (via dl[, 1L] in to_install) so R CMD check --as-cran can run, but we
## do not chase Suggests recursively.
cat("\n--- Resolving and installing transitive dependencies ---\n")
with_suggests <- identical(Sys.getenv("WITH_SUGGESTS", "false"), "true")
dep_types <- if (with_suggests) {
    c("Depends", "Imports", "LinkingTo", "Suggests")
} else {
    c("Depends", "Imports", "LinkingTo")
}
cat(sprintf("Recursive dependency types: %s\n", paste(dep_types, collapse = ", ")))
all_deps <- tools::package_dependencies(
    dl[, 1L], db = ap,
    which     = dep_types,
    recursive = TRUE)
to_install <- sort(unique(c(dl[, 1L],
                             unlist(all_deps, use.names = FALSE))))

## Drop base/recommended packages already present in every R installation,
## and lme4 itself (installed separately into versioned library dirs below)
already <- rownames(installed.packages(
    priority = c("base", "recommended")))
to_install <- setdiff(to_install, c(already, "lme4"))
cat(sprintf("Installing %d packages ...\n", length(to_install)))

## configure.vars for packages with non-standard build systems
## (mirrors the configure.vars block in checkReverse.R)
configure.vars <- list(
    arrow = c("LIBARROW_BINARY=false",
              "ARROW_R_DEV=true",
              "ARROW_DEPENDENCY_SOURCE=BUNDLED"))

## Try a fast batch install first.  If a single package (e.g. r-cran-micemd)
## has a broken post-install script it will leave dpkg in an interrupted state
## and fail the whole call.  In that case, clean up the dpkg state and fall
## back to per-package installs so individual failures can be skipped.
##
## Note: bspm intercepts install.packages() and serializes installs through
## apt, so Ncpus only affects within-package compilation, not across-package
## parallelism.  r2u covers most of CRAN as binaries but only a subset of
## Bioconductor, so source compilation can be a bottleneck when Bioc packages
## are included.  If build times are unacceptable, consider doing a binary-only
## pass first (bspm active), then calling bspm::disable() before a second pass
## with install.packages(..., Ncpus = NCPUS) for remaining source-only packages.
batch_ok <- tryCatch({
    install.packages(to_install,
                     configure.vars = configure.vars,
                     dependencies   = FALSE,
                     Ncpus          = NCPUS)
    TRUE
}, error = function(e) {
    message("Batch install failed: ", conditionMessage(e))
    FALSE
})

if (!batch_ok) {
    message("Fixing broken dpkg state and retrying per-package ...")
    system2("dpkg", c("--configure", "-a"))
    system2("apt-get", c("install", "-f", "-y"))

    install_failed <- character(0)
    still_needed <- setdiff(to_install, rownames(installed.packages()))
    cat(sprintf("Retrying %d packages individually ...\n", length(still_needed)))
    for (pkg in still_needed) {
        tryCatch(
            install.packages(pkg, dependencies = FALSE, Ncpus = NCPUS),
            error = function(e) {
                message(sprintf("  SKIP %s: %s", pkg, conditionMessage(e)))
                install_failed <<- c(install_failed, pkg)
            }
        )
    }
    if (length(install_failed)) {
        warning(sprintf("%d package(s) could not be installed: %s",
                        length(install_failed),
                        paste(install_failed, collapse = ", ")))
        writeLines(install_failed,
                   file.path(REVDEP_DIR, "install_failures.txt"))
    }
}

## ---- 5 & 6. Install old and new lme4 into separate library dirs ----------
## At check time, check_one.R prepends the appropriate library so that the
## right lme4 version is found first, without rebuilding the whole image.
cat("\n--- Installing lme4 versions into versioned library dirs ---\n")
for (info in list(c(lme4_old, LIB_OLD), c(lme4_new, LIB_NEW))) {
    tgz <- info[1]; lib <- info[2]
    cat(sprintf("Installing %s into %s ...\n", basename(tgz), lib))
    install.packages(tgz, lib = lib, repos = NULL,
                     type = "source", Ncpus = NCPUS)
}

## ---- Summary -------------------------------------------------------------
n_check <- length(readLines(PKG_LIST_FILE))
summary_lines <- c(
    sprintf("old lme4       : %s  -> %s", basename(lme4_old), LIB_OLD),
    sprintf("new lme4       : %s  -> %s", basename(lme4_new), LIB_NEW),
    sprintf("Revdeps found  : %d", length(rdeps)),
    sprintf("Tarballs saved : %d  (see %s)", n_check, PKG_LIST_FILE),
    sprintf("Install failures: %s",
            if (length(failed)) paste(failed, collapse = ", ") else "none"),
    sprintf("R library      : %s", .libPaths()[[1L]])
)
writeLines(summary_lines, file.path(REVDEP_DIR, "setup_summary.txt"))
cat("\n=== Setup complete ===\n")
writeLines(summary_lines)
