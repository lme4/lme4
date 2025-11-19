checkReverse <-
function (args) {

    args.prefix <- sub("^(--.*?=)(.*)$", "\\1", args)
    args.suffix <- sub("^(--.*?=)(.*)$", "\\2", args)

    filename <- args[i.filename <- grep("^[^-].*[.]tar[.]gz$", args)]
    stopifnot(length(filename) == 1L, file.exists(filename))

    libclean <- FALSE
    preclean <- TRUE
       clean <- FALSE
    libpaths <- .Library
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
            "--library=" =
                libpaths <- c(strsplit(args.suffix[[i]], ":")[[1L]],
                              libpaths),
            "--jobs=" =
                   Ncpus <- as.integer(args.suffix[[i]]),
            if (i != i.filename)
            stop(gettextf("invalid command line argument '%s'",
                          args[[i]]),
                 domain = NA))

    .lp <- .libPaths()
    on.exit(.libPaths(.lp))
    .libPaths(libpaths)

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

    repos <- utils:::.expand_BioC_repository_URLs(
        c("CRAN"     = "https://cloud.r-project.org",
          "BioCsoft" = "%bm/packages/%v/bioc",
          "BioCann"  = "%bm/packages/%v/data/annotation",
          "BioCexp"  = "%bm/packages/%v/data/experiment",
          "INLA"     = "https://inla.r-inla-download.org/R/stable",
          "CmdStan"  = "https://stan-dev.r-universe.dev"))
    .op <- options(repos = repos)
    on.exit(options(.op), add = TRUE)

    reverse <- list(repos = repos,
                    which = c("Depends", "Imports", "LinkingTo"),
                    recursive = FALSE)

    ans <- tools::check_packages_in_dir(dirname,
                                        reverse = reverse,
                                        Ncpus = Ncpus,
                                        clean = clean)
    ans
}

args <- commandArgs(trailingOnly = TRUE)
ch <- checkReverse(args)


## 'tkrplot', 'arrow'
