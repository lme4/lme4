## Usage in a shell:
##
## $ R -f reverse.R --args [option1] ... [optionN] /path/to/tarball
##
## Usage in R:
##
## > checkReverse(c("[option1]", ..., "[optionN]", "/path/to/tarball"))
##
## Example:
##
## $ LME4_OLD=lme4_1.1-37.tar.gz
## $ LME4_NEW=lme4_1.1-38.tar.gz
## $ R --vanilla -f reverse.R --args --jobs=4 ${LME4_OLD}
## $ R --vanilla -f reverse.R --args --jobs=4 ${LME4_NEW}
## $ R -e "tools::check_packages_in_dir_changes(\"${LME4_NEW}.reverse\", \"${LME4_OLD}.reverse\", TRUE, TRUE)[\"<\", ]" # changes to worse
##
checkReverse <-
function (args) {
    ## Wishlist for KH/tools::check_packages_in_dir:
    ## * support no-check/only-install mode
    ## * support building library tree on top of package being checked

    stopifnot(is.character(args))
    args.prefix <- sub("^(--.*?=)(.*)$", "\\1", args)
    args.suffix <- sub("^(--.*?=)(.*)$", "\\2", args)

    filename <- args[i.filename <- grep("^[^-].*[.]tar[.]gz$", args)]
    stopifnot(length(filename) == 1L, file.exists(filename))

    libclean <- FALSE
    preclean <- TRUE
       clean <- FALSE
       check <- TRUE
    libpaths <- character(0L)
       Ncpus <- 1L

    for (i in seq_along(args))
    switch (args.prefix[[i]],
            "--libclean" =
                libclean <- TRUE,
            "--no-libclean" =
                libclean <- FALSE,
            "--preclean" =
                preclean <- TRUE,
            "--no-preclean" =
                preclean <- FALSE,
            "--clean" =
                   clean <- TRUE,
            "--no-clean" =
                   clean <- FALSE,
            "--check" =
                   check <- TRUE,
            "--no-check" =
                   check <- FALSE,
            "--library=" =
                libpaths <- strsplit(args.suffix[[i]], ":")[[1L]],
            "--jobs=" =
                   Ncpus <- as.integer(args.suffix[[i]]),
            if (i != i.filename)
            stop(gettextf("invalid command line argument '%s'",
                          args[[i]]),
                 domain = NA))

    dirname <- paste0(filename, ".reverse")
    if (!dir.exists(dirname))
        dir.create(dirname)
    else {
        if (libclean)
            unlink(file.path(dirname, "Library"),
                   recursive = TRUE)
        if (preclean)
            unlink(file.path(dirname, c("Outputs",
                                        "*.Rcheck",
                                        "*.tar.gz",
                                        "PACKAGES",
                                        "PACKAGES.gz",
                                        "PACKAGES.rds",
                                        "timings.tab")),
                   recursive = TRUE)
    }
    file.copy(filename, dirname, overwrite = TRUE)

    .lp <- .libPaths()
    on.exit(.libPaths(.lp), add = TRUE)
    .libPaths(libpaths)

    repos <- utils:::.expand_BioC_repository_URLs(
        c("CRAN"     = "https://cloud.r-project.org",
          "BioCsoft" = "%bm/packages/%v/bioc",
          "BioCann"  = "%bm/packages/%v/data/annotation",
          "BioCexp"  = "%bm/packages/%v/data/experiment",
          "INLA"     = "https://inla.r-inla-download.org/R/stable",
          "CmdStan"  = "https://stan-dev.r-universe.dev"))
    .op <- options(repos = repos, useFancyQuotes = FALSE)
    on.exit(options(.op), add = TRUE)

    cpid <- tools::check_packages_in_dir
    if (!check) {
        ## Insert early return after call to utils::install.packages;
        ## the return value is a character vector storing the names of
        ## packages whose installation failed.
        stopifnot(getRversion() >= "4.3", getRversion() < "4.7")
        body(cpid)[[43L]] <-
        substitute(env = list(.__ORIG__. = body(cpid)[[43L]]), {
            outfiles <- Sys.glob(file.path(outdir, "install_*.out"))
            outnames <- sub("^install_(.*)[.]out$", "\\1", basename(outfiles))
            grrl <- function (out) grep("^[*] removing", readLines(out))
            outfails <- lengths(lapply(outfiles, grrl)) > 0L
            return(outnames[outfails])
            .__ORIG__.
        })
    }
    cpid(dirname, reverse = list(), Ncpus = Ncpus, clean = clean)
}

args <- commandArgs(trailingOnly = TRUE)
ch <- checkReverse(args)
