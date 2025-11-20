## USAGE IN A SHELL
##
##     $ R -f checkChanges.R [option1] ... [optionN] --old=<directory> --new=<directory>
##
## USAGE IN R
##
##     > checkChanges(c("[option1]", ..., "[optionN]", "--old=<directory>", "--new=<directory>"))
##
## EXAMPLE
##
##     $ LME4_OLD=lme4_1.1-37.tar.gz
##     $ LME4_NEW=lme4_1.1-38.tar.gz
##     $ R --vanilla -f checkReverse.R --args --jobs=4 ${LME4_OLD}
##     $ R --vanilla -f checkReverse.R --args --jobs=4 ${LME4_NEW}
##     $ R --vanilla -f checkChanges.R --args --old=${LME4_OLD}.reverse --new=${LME4_NEW}.reverse
##
checkChanges <-
function (args) {
    stopifnot(is.character(args))
    args.prefix <- sub("^(--.*?=)(.*)$", "\\1", args)
    args.suffix <- sub("^(--.*?=)(.*)$", "\\2", args)

    preclean <- TRUE
    olddir <- newdir <- outdir <- NULL

    for (i in seq_along(args))
    switch (args.prefix[[i]],
            "--preclean" =
                preclean <- TRUE,
            "--no-preclean" =
                preclean <- FALSE,
            "--old=" =
                  olddir <- args.suffix[[i]],
            "--new=" =
                  newdir <- args.suffix[[i]],
            "--output=" =
                  outdir <- args.suffix[[i]],
            stop(gettextf("invalid command line option '%s'",
                          args[[i]]),
                 domain = NA))

    stopifnot(!is.null(olddir), !is.null(newdir))
    if (is.null(outdir))
        outdir <- sprintf("CHANGES_old=%s_new=%s", olddir, newdir)
    if (!dir.exists(outdir))
        dir.create(outdir)
    else {
        if (preclean)
            unlink(file.path(outdir, c("CHANGES",
                                       "CHANGES.rds",
                                       "*_00install.out",
                                       "*_00check.log")))
    }

    changes <- tools::check_packages_in_dir_changes(newdir, olddir, TRUE, TRUE)
    changes <- changes["<", ] # to worse
    package <- unique(changes[["Package"]])

    capture.output(print(changes), file = file.path(outdir, "CHANGES"),
                   split = TRUE)
    saveRDS(changes, file = file.path(outdir, "CHANGES.rds"))
    for (zz in c("00install.out", "00check.log"))
    file.copy(file.path(newdir,
                        sprintf("rdepends_%s.Rcheck", package),
                        zz),
              file.path(outdir,
                        sprintf("%s_%s", package, zz)),
              overwrite = TRUE)

    changes
}

args <- commandArgs(trailingOnly = TRUE)
ch <- checkChanges(args)
