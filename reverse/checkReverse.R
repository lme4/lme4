## USAGE IN A SHELL
##
##     $ R -f checkReverse.R --args [option1] ... [optionN] <archive>
##
## USAGE IN R
##
##     > checkReverse(c("[option1]", ..., "[optionN]", "<archive>"))
##
## EXAMPLE
##
##     Run without the check until you no longer encounter spurious
##     install failures.  Packages may require manual intervention.
##     Once you are satisfied, run with the check.
##
##     $ LME4_ARCHIVE=lme4_1.1-37.tar.gz
##     $ R --vanilla -f checkReverse.R --args --jobs=4 --no-check ${LME4_ARCHIVE}
##     $ R --vanilla -f checkReverse.R --args --jobs=4            ${LME4_ARCHIVE}
##
##     If you want to make use of an existing library tree for a given
##     run, then use --library=, giving a colon-delimited search path.
##
checkReverse <-
function (args) {
    ## Wishlist for KH/tools::check_packages_in_dir:
    ## * support no-check/only-install mode
    ## * support building library tree on top of package being checked

    stopifnot(is.character(args))
    args.prefix <- sub("^(--.*?=)(.*)$", "\\1", args)
    args.suffix <- sub("^(--.*?=)(.*)$", "\\2", args)

    tarfile <- args[i.tarfile <- grep("^[^-].*[.]tar[.]gz$", args)]
    stopifnot(length(tarfile) == 1L)

    libclean <- FALSE
    preclean <- TRUE
       clean <- FALSE
       check <- TRUE
      outdir <- NULL
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
            "--output=" =
                  outdir <- args.suffix[[i]],
            "--library=" =
                libpaths <- strsplit(args.suffix[[i]], ":")[[1L]],
            "--jobs=" =
                   Ncpus <- as.integer(args.suffix[[i]]),
            if (i != i.tarfile)
            stop(gettextf("invalid command line option '%s'",
                          args[[i]]),
                 domain = NA))

    if (is.null(outdir))
        outdir <- sprintf("%s.reverse", tarfile)
    if (!dir.exists(outdir))
        dir.create(outdir)
    else {
        if (libclean)
            unlink(file.path(outdir, "Library"),
                   recursive = TRUE)
        if (preclean)
            unlink(file.path(outdir, c("Outputs",
                                       "*.Rcheck",
                                       "*.tar.gz",
                                       "PACKAGES",
                                       "PACKAGES.gz",
                                       "PACKAGES.rds",
                                       "timings.tab")),
                   recursive = TRUE)
    }
    stopifnot(file.copy(tarfile, outdir, overwrite = TRUE))

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
    .op <- options(repos = repos,
                   install.packages.compile.from.source = TRUE,
                   useFancyQuotes = FALSE)
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
    configure.args <-
        list()
    configure.vars <-
        list(arrow = c("LIBARROW_BINARY=false",
                       "ARROW_R_DEV=true",
                       "ARROW_DEPENDENCY_SOURCE=BUNDLED",
                       if (Sys.info()[["user"]] == "mikael")
                           "PKG_CONFIG=\"pkg-config --static\"",
                       NULL))
    out <- cpid(outdir, reverse = list(), Ncpus = Ncpus, clean = clean,
                install_args = list(configure.args = configure.args,
                                    configure.vars = configure.vars))
    saveRDS(out, file = file.path(outdir, "REVERSE.rds"))
    out
}

args <- commandArgs(trailingOnly = TRUE)
ch <- checkReverse(args)
