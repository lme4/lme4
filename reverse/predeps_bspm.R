## predeps_bspm.R -- pre-install, via r2u/bspm binaries, the dependency
## closure of a set of package tarballs, so that a later
## tools::check_packages_in_dir() run (which always installs with
## type="source") finds them already satisfied and has little/nothing left
## to compile.
##
## Usage:
##   R --vanilla -f predeps_bspm.R --args pkglist.txt
##
## Installs into the default (site) library, i.e. wherever bspm/apt puts
## r2u binaries -- NOT into a per-check Library/ subdirectory, so this is
## meant to be run once up front, shared across check runs.

args <- commandArgs(trailingOnly = TRUE)
pkglistfile <- args[[1L]]
pkgs <- readLines(pkglistfile)
pkgs <- trimws(pkgs)
pkgs <- pkgs[nzchar(pkgs) & !grepl("^#", pkgs)]
cat(sprintf("Primary packages (%d): %s\n", length(pkgs), paste(pkgs, collapse = ", ")))

stopifnot(requireNamespace("bspm", quietly = TRUE))
suppressMessages(bspm::enable())
options(repos = c(CRAN = "https://cloud.r-project.org"),
        useFancyQuotes = FALSE)

NCPUS <- max(1L, parallel::detectCores() - 2L)

cat("--- Resolving dependency closure (Depends/Imports/LinkingTo, recursive) ---\n")
ap <- available.packages()
hard <- tools::package_dependencies(pkgs, db = ap,
                                     which = c("Depends", "Imports", "LinkingTo"),
                                     recursive = TRUE)
## One level of Suggests for the primary packages themselves (needed to run
## their example/test/vignette suites under --as-cran), not chased recursively.
suggests1 <- tools::package_dependencies(pkgs, db = ap, which = "Suggests",
                                          recursive = FALSE)

to_install <- sort(unique(c(pkgs,
                            unlist(hard, use.names = FALSE),
                            unlist(suggests1, use.names = FALSE))))

already <- rownames(installed.packages(priority = c("base", "recommended")))
to_install <- setdiff(to_install, c(already, "lme4"))
cat(sprintf("Installing %d package(s) via bspm/r2u (type=\"both\") ...\n", length(to_install)))
print(to_install)

install.packages(to_install, type = "both", dependencies = FALSE, Ncpus = NCPUS)

still_missing <- setdiff(to_install, rownames(installed.packages()))
cat(sprintf("\n%d package(s) still missing after binary pass (will fall back to source in check step):\n",
            length(still_missing)))
print(still_missing)
