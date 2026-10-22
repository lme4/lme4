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
##
## With REFRESH=true (Dockerfile.refresh, 'build.sh --refresh BASE'), update
## an image previously built by this script instead of starting from
## scratch: step 3 only downloads tarballs that are new or have a new version
## (and deletes stale ones), and step 4 first upgrades what is already
## installed (apt for r2u binaries, source reinstalls for packages in the
## source-built library) before installing whatever is newly needed.

options(timeout = 120)   # default 60s can time out on slow connections

REVDEP_DIR    <- "/opt/revdep"
TARBALL_DIR   <- file.path(REVDEP_DIR, "tarballs")
PKG_LIST_FILE <- file.path(REVDEP_DIR, "pkgs_to_check.txt")
LIB_OLD       <- file.path(REVDEP_DIR, "Library_old")
LIB_NEW       <- file.path(REVDEP_DIR, "Library_new")
NCPUS         <- max(1L, parallel::detectCores() - 1L)
refresh       <- identical(Sys.getenv("REFRESH", "false"), "true")
## where pass 2 (source compilation) installs; r2u binaries go to
## /usr/lib/R/site-library
SRC_LIB       <- .libPaths()[[1L]]

## stale from a base image's build
unlink(file.path(REVDEP_DIR, "install_failures.txt"))

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
##
## When refreshing, tarballs whose file name (PKG_VERSION.tar.gz) matches the
## current repo version are kept; the rest of the directory (old versions,
## packages that are no longer reverse dependencies) is deleted.
cat("\n--- Downloading source tarballs ---\n")
want_tgz <- paste0(rdeps, "_", ap[rdeps, "Version"], ".tar.gz")
names(want_tgz) <- rdeps
n_removed <- 0L
if (refresh) {
    have_tgz <- list.files(TARBALL_DIR, pattern = "\\.tar\\.gz$")
    stale <- setdiff(have_tgz, want_tgz)
    unlink(file.path(TARBALL_DIR, stale))
    n_removed <- length(stale)
    to_download <- rdeps[!want_tgz %in% have_tgz]
    cat(sprintf("Refresh: keeping %d tarball(s), removed %d stale, downloading %d\n",
                length(rdeps) - length(to_download), n_removed,
                length(to_download)))
} else {
    to_download <- rdeps
}
dl <- download.packages(to_download, destdir = TARBALL_DIR,
                        repos = repos, type = "source")
failed <- setdiff(to_download, dl[, 1L])
if (length(failed))
    warning(sprintf("Could not download %d package(s): %s",
                    length(failed), paste(failed, collapse = ", ")))
cat(sprintf("Downloaded %d / %d tarballs to %s\n",
            nrow(dl), length(to_download), TARBALL_DIR))
## all revdeps whose tarball is present, whether downloaded now or kept
rdeps_ok <- setdiff(rdeps, failed)
writeLines(file.path(TARBALL_DIR, want_tgz[rdeps_ok]), PKG_LIST_FILE)

## ---- 4. Install all transitive dependencies (excluding lme4) -------------
## Resolve the full dependency closure of the downloaded rev-deps so that
## R CMD check can run offline on Compute Canada.
## Two-pass: r2u/bspm installs binaries first, then source compilation in parallel.
##
## Recursive resolution uses only hard deps (Depends/Imports/LinkingTo) to
## avoid pulling in the enormous transitive closure of all Suggests.
## Suggests of the rev-dep packages themselves are included one level deep
## (via rdeps_ok in to_install) so R CMD check --as-cran can run, but we
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
    rdeps_ok, db = ap,
    which     = dep_types,
    recursive = TRUE)
to_install <- sort(unique(c(rdeps_ok,
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

## ---- 4a. (refresh only) upgrade what is already installed -----------------
## r2u binaries are ordinary .debs, so apt fetches only those that changed.
## r-base-core is held: a new R minor version would invalidate every
## installed package, and then a full rebuild is the right tool.
## Source-built packages (mostly Bioconductor) are updated in place in
## SRC_LIB, with bspm disabled so that a newer r2u binary installed into
## /usr/lib/R/site-library can't end up shadowed by the outdated source
## build. Source-built packages that were not themselves updated but have
## a LinkingTo dependency that was (Rcpp, RcppEigen, TMB, ...) were compiled
## against the old headers, so they are reinstalled too.
## installed packages, keeping only the copy that loads first when a package
## is in more than one library
installed_first <- function() {
    ip <- installed.packages()
    ip[!duplicated(ip[, "Package"]), , drop = FALSE]
}
changed_since <- function(ver) {
    now <- installed_first()[, "Version"]
    names(ver)[ver != now[names(ver)]]
}
to_install_now <- to_install
n_apt <- n_src <- n_relink <- 0L
if (refresh) {
    cat("\n--- Refresh: upgrading installed packages ---\n")
    ver_before <- installed_first()[, "Version"]
    system2("apt-mark", c("hold", "r-base-core"))
    stopifnot(system2("apt-get", "update") == 0L,
              system2("apt-get", c("upgrade", "-y", "--with-new-pkgs",
                                   "--no-install-recommends")) == 0L)
    n_apt <- length(changed_since(ver_before))
    cat(sprintf("apt upgraded %d R package(s)\n", n_apt))

    bspm::disable()
    op <- old.packages(lib.loc = SRC_LIB, repos = repos)
    if (!is.null(op)) {
        n_src <- nrow(op)
        cat(sprintf("Updating %d source-built package(s): %s\n",
                    n_src, paste(rownames(op), collapse = ", ")))
        install.packages(rownames(op), lib = SRC_LIB,
                         configure.vars = configure.vars,
                         dependencies = FALSE, Ncpus = NCPUS)
    }
    changed <- changed_since(ver_before)
    src_pkgs <- setdiff(rownames(installed.packages(lib.loc = SRC_LIB)),
                        changed)
    links <- tools::package_dependencies(src_pkgs, db = installed_first(),
                                         which = "LinkingTo")
    relink <- names(Filter(function(d) any(d %in% changed), links))
    n_relink <- length(relink)
    if (n_relink) {
        cat(sprintf("Reinstalling %d package(s) with upgraded LinkingTo deps: %s\n",
                    n_relink, paste(relink, collapse = ", ")))
        install.packages(relink, lib = SRC_LIB,
                         configure.vars = configure.vars,
                         dependencies = FALSE, Ncpus = NCPUS)
    }
    bspm::enable()
    to_install_now <- setdiff(to_install, rownames(installed.packages()))
    cat(sprintf("%d newly needed package(s) to install\n",
                length(to_install_now)))
}

## Two-pass install strategy:
##
## Pass 1 (bspm active): install.packages() is intercepted by bspm, which
##   installs pre-built binary .deb packages via apt where available (fast).
##   bspm silently skips packages with no pre-built binary.  If a single
##   package has a broken post-install script and leaves dpkg in a bad state,
##   we detect and recover before proceeding.
##   Observed timing: ~25 min (standard build, WITH_BIOC=true, WITH_SUGGESTS=false).
##
## Pass 2 (bspm disabled): find packages still missing after pass 1 and
##   compile them from source using R's own parallel installer (Ncpus
##   subprocesses).  r2u covers most of CRAN as binaries but only a subset of
##   Bioconductor, so this pass handles source-only packages in parallel
##   rather than letting bspm serialize them one at a time through apt.
##   Observed timing: ~8 min (standard build, WITH_BIOC=true, WITH_SUGGESTS=false).
cat("\n--- Pass 1: binary installs via bspm/apt ---\n")
batch_ok <- length(to_install_now) == 0L || tryCatch({
    install.packages(to_install_now,
                     configure.vars = configure.vars,
                     dependencies   = FALSE,
                     Ncpus          = NCPUS)
    TRUE
}, error = function(e) {
    message("Pass 1 batch install failed: ", conditionMessage(e))
    FALSE
})

if (!batch_ok) {
    message("Fixing broken dpkg state and retrying pass 1 ...")
    system2("dpkg", c("--configure", "-a"))
    system2("apt-get", c("install", "-f", "-y"))
    tryCatch(
        install.packages(to_install_now, configure.vars = configure.vars,
                         dependencies = FALSE, Ncpus = NCPUS),
        error = function(e) message("Pass 1 retry also failed: ", conditionMessage(e))
    )
}

## Pass 2: parallel source compilation for packages bspm could not supply.
## install.packages() only warns (does not error) when an individual package
## fails to install, so failures are detected by checking what is actually
## installed afterwards, not by catching errors.
still_needed <- setdiff(to_install, rownames(installed.packages()))
if (length(still_needed) > 0L) {
    cat(sprintf("\n--- Pass 2: source compilation for %d package(s) (Ncpus=%d) ---\n",
                length(still_needed), NCPUS))
    bspm::disable()
    tryCatch(
        install.packages(still_needed,
                         configure.vars = configure.vars,
                         dependencies   = FALSE,
                         Ncpus          = NCPUS),
        error = function(e) message("Pass 2 parallel install failed: ", conditionMessage(e))
    )
    still_needed2 <- setdiff(still_needed, rownames(installed.packages()))
    if (length(still_needed2) > 0L) {
        cat(sprintf("Retrying %d package(s) individually ...\n", length(still_needed2)))
        for (pkg in still_needed2) {
            tryCatch(
                install.packages(pkg, configure.vars = configure.vars,
                                 dependencies = FALSE, Ncpus = NCPUS),
                error = function(e)
                    message(sprintf("  SKIP %s: %s", pkg, conditionMessage(e)))
            )
        }
    }
    bspm::enable()
} else {
    cat("Pass 2: no source-only packages remaining.\n")
}

install_failed <- setdiff(to_install, rownames(installed.packages()))
if (length(install_failed) > 0L) {
    warning(sprintf("%d package(s) could not be installed: %s",
                    length(install_failed),
                    paste(install_failed, collapse = ", ")))
    writeLines(install_failed,
               file.path(REVDEP_DIR, "install_failures.txt"))
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
    sprintf("Download failures: %s",
            if (length(failed)) paste(failed, collapse = ", ") else "none"),
    sprintf("Install failures: %s",
            if (length(install_failed)) paste(install_failed, collapse = ", ") else "none"),
    sprintf("R library      : %s", .libPaths()[[1L]])
)
if (refresh)
    summary_lines <- c(summary_lines,
        sprintf("refreshed      : %s; tarballs downloaded %d, removed %d; R pkgs upgraded via apt %d, from source %d, relinked %d",
                format(Sys.Date()), nrow(dl), n_removed, n_apt, n_src, n_relink))
writeLines(summary_lines, file.path(REVDEP_DIR, "setup_summary.txt"))
cat("\n=== Setup complete ===\n")
writeLines(summary_lines)
