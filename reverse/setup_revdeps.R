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
## Use BiocManager for a consistent CRAN + Bioconductor URL set.
## (mirrors utils:::.expand_BioC_repository_URLs used in checkReverse.R but
##  without depending on an internal function)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
repos <- BiocManager::repositories()
options(repos = repos,
        Ncpus = NCPUS,
        install.packages.compile.from.source = "always",
        useFancyQuotes = FALSE)
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
all_deps <- tools::package_dependencies(
    dl[, 1L], db = ap,
    which     = c("Depends", "Imports", "LinkingTo"),
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

install.packages(to_install,
                 configure.vars = configure.vars,
                 dependencies   = FALSE,   # closure already in to_install
                 Ncpus          = NCPUS)

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
