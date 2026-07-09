## check_local_subset.R -- download a named list of CRAN packages and run
## R CMD check on them using whatever lme4 is currently installed on
## .libPaths() (no reverse-dependency discovery, no separate lme4 library).
##
## Usage:
##   R --vanilla -f check_local_subset.R --args --jobs=8 --output=DIR pkglist.txt
##
## pkglist.txt: one package name per line (blank lines/lines starting with
## '#' ignored).

args <- commandArgs(trailingOnly = TRUE)
args.prefix <- sub("^(--.*?=)(.*)$", "\\1", args)
args.suffix <- sub("^(--.*?=)(.*)$", "\\2", args)

pkglistfile <- args[grep("^[^-]", args)]
stopifnot(length(pkglistfile) == 1L)

Ncpus  <- 1L
outdir <- NULL
for (i in seq_along(args)) {
    if (args[[i]] == pkglistfile) next
    switch(args.prefix[[i]],
           "--jobs="   = Ncpus  <- as.integer(args.suffix[[i]]),
           "--output=" = outdir <- args.suffix[[i]],
           stop(gettextf("invalid command line option '%s'", args[[i]]), domain = NA))
}
if (is.null(outdir)) outdir <- "local_check"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pkgs <- readLines(pkglistfile)
pkgs <- trimws(pkgs)
pkgs <- pkgs[nzchar(pkgs) & !grepl("^#", pkgs)]
cat(sprintf("Packages to check (%d): %s\n", length(pkgs), paste(pkgs, collapse = ", ")))

cat(sprintf("lme4 installed: %s (from %s)\n",
            as.character(packageVersion("lme4")),
            dirname(find.package("lme4"))))

options(repos = c(CRAN = "https://cloud.r-project.org"),
        install.packages.compile.from.source = TRUE,
        useFancyQuotes = FALSE)

## Force r2u/bspm binary installs even when invoked with --vanilla (which
## skips /etc/R/Rprofile.site, where bspm::enable() normally happens).
if (requireNamespace("bspm", quietly = TRUE)) {
    suppressMessages(bspm::enable())
    cat("bspm enabled: using r2u binary packages where available\n")
} else {
    cat("bspm not installed: falling back to source installs\n")
}

## ---- Download source tarballs into outdir --------------------------------
existing <- Sys.glob(file.path(outdir, "*.tar.gz"))
have <- sub("_.*$", "", basename(existing))
need <- setdiff(pkgs, have)
if (length(need)) {
    cat(sprintf("Downloading %d tarball(s): %s\n", length(need), paste(need, collapse = ", ")))
    dl <- download.packages(need, destdir = outdir, type = "source")
    failed <- setdiff(need, dl[, 1L])
    if (length(failed))
        warning(sprintf("Could not download: %s", paste(failed, collapse = ", ")))
} else {
    cat("All tarballs already present.\n")
}

## ---- Run check_packages_in_dir on just these tarballs --------------------
## pfiles must be given as basenames (relative to dir); check_packages_in_dir
## internally does file.path(dir, pfiles) itself.
pfiles <- basename(Sys.glob(file.path(outdir, "*.tar.gz")))
pfiles <- pfiles[sub("_.*$", "", pfiles) %in% pkgs]
cat(sprintf("Checking %d tarball(s):\n", length(pfiles)))
print(pfiles)

out <- tools::check_packages_in_dir(
    outdir,
    pfiles = pfiles,
    check_args = "--as-cran",
    check_env = c("_R_CHECK_FORCE_SUGGESTS_=false"),
    Ncpus = Ncpus,
    clean = FALSE)

saveRDS(out, file = file.path(outdir, "LOCAL_CHECK.rds"))
tryCatch(print(summary(out)),
         error = function(e)
             cat(sprintf("(summary() failed: %s -- see .Rcheck dirs directly)\n",
                          conditionMessage(e))))
