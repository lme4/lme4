#!/usr/bin/env Rscript
##
## setup_revdeps.R -- internet-connected setup phase, run during 'docker build'
##
## Adapted from checkReverse.R (the install-only / --no-check phase).
##
## Steps:
##   1. Configure CRAN + Bioconductor repositories
##   2. Discover all direct reverse dependencies of lme4
##   3. Download their source tarballs into /opt/revdep/tarballs/
##   4. Install all transitive dependencies into the container R library
##   5. Install dev lme4 from its tarball (overwrites any CRAN binary)

REVDEP_DIR    <- "/opt/revdep"
TARBALL_DIR   <- file.path(REVDEP_DIR, "tarballs")
PKG_LIST_FILE <- file.path(REVDEP_DIR, "pkgs_to_check.txt")
NCPUS         <- max(1L, parallel::detectCores())

dir.create(TARBALL_DIR, recursive = TRUE, showWarnings = FALSE)

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

## ---- 4. Install all transitive dependencies ------------------------------
## Resolve the full dependency closure of the downloaded rev-deps so that
## R CMD check can run offline on Compute Canada.
## r2u installs binary apt packages where available (fast); falls back to source.
cat("\n--- Resolving and installing transitive dependencies ---\n")
all_deps <- tools::package_dependencies(
    dl[, 1L], db = ap,
    which     = c("Depends", "Imports", "LinkingTo", "Suggests"),
    recursive = TRUE)
to_install <- sort(unique(c(dl[, 1L],
                             unlist(all_deps, use.names = FALSE))))

## Drop base/recommended packages already present in every R installation
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

## ---- 5. Install dev lme4 from tarball ------------------------------------
cat("\n--- Installing dev lme4 from source tarball ---\n")
lme4_tgz <- Sys.glob(file.path(REVDEP_DIR, "lme4_*.tar.gz"))
stopifnot("exactly one lme4 tarball expected" = length(lme4_tgz) == 1L)
install.packages(lme4_tgz, repos = NULL, type = "source", Ncpus = NCPUS)
cat(sprintf("Installed lme4 from %s\n", basename(lme4_tgz)))

## ---- Summary -------------------------------------------------------------
n_check <- length(readLines(PKG_LIST_FILE))
summary_lines <- c(
    sprintf("lme4 source    : %s", basename(lme4_tgz)),
    sprintf("Revdeps found  : %d", length(rdeps)),
    sprintf("Tarballs saved : %d  (see %s)", n_check, PKG_LIST_FILE),
    sprintf("Install failures: %s",
            if (length(failed)) paste(failed, collapse = ", ") else "none"),
    sprintf("R library      : %s", .libPaths()[[1L]])
)
writeLines(summary_lines, file.path(REVDEP_DIR, "setup_summary.txt"))
cat("\n=== Setup complete ===\n")
writeLines(summary_lines)
