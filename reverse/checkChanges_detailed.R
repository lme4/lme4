## checkChanges_detailed.R -- like checkChanges.R, but compares the *messages*
## of each check (ERROR/WARNING/NOTE output), not just its status.
##
## checkChanges.R keeps only checks whose status got worse (e.g. OK -> ERROR),
## so a check that already failed under the old lme4 for an unrelated reason
## (no network on compute nodes, a Suggests package missing from the image,
## ...) hides any new lme4-related failure in the same check: ERROR -> ERROR
## counts as "no change". This script classifies every check that has an
## issue under either version as
##
##   worse   : status got worse          (what checkChanges.R reports)
##   changed : same status, but the (normalized) messages differ
##   better  : status got better         (only reported with --include-better)
##
## and, for each, lists the message lines that appear only under the old or
## only under the new version.
##
## Before comparing, outputs are normalized to drop run-to-run noise:
## fancy quotes, timings such as [35s/39s], the "Examples with CPU ... > 5s"
## timing tables, R temp-file names, and the .Rcheck directory prefix.
## Other numbers are kept, so changed numerical test results still show up.
##
## USAGE IN A SHELL
##
##     $ R -f checkChanges_detailed.R --args [option1] ... [optionN] --old=<directory> --new=<directory>
##
## OPTIONS
##
##     --old=DIR, --new=DIR : directories of rdepends_PKG.Rcheck outputs
##     --output=DIR         : where to write the results (default:
##                            CHANGES_DETAILED_old=<old>_new=<new>)
##     --include-better     : also report checks whose status improved
##     --full               : print full old/new outputs, not just the
##                            lines that differ
##     --max-lines=N        : show at most the last N differing lines per
##                            side of each check (default 30; 0 = all).
##                            Tail rather than head, since errors and test
##                            summaries come at the end of the output.
##     --[no-]preclean      : remove previous results from DIR (default yes)
##
## Writes, in the output directory:
##     CHANGES_detailed      : text report
##     CHANGES_detailed.rds  : data frame (Package, Check, Type, Old, New,
##                             Removed, Added, OldOutput, NewOutput)
##     PKG_00install.out, PKG_00check.log : new-version logs for each
##                             package with a reported change
##
## EXAMPLE (both-mode results from slurm_submit.sh)
##
##     $ R --vanilla -f checkChanges_detailed.R --args \
##           --old=results/old --new=results/new --output=CHANGES_detailed_run1

issue_levels <- c("", "OK", "NOTE", "WARNING", "ERROR", "FAILURE")

## rank of a check status; statuses outside issue_levels (e.g. INFO) and
## missing checks rank as ""
status_rank <- function(s) {
    s[is.na(s) | is.na(match(s, issue_levels))] <- ""
    match(s, issue_levels)
}

normalize_output <- function(x) {
    x[is.na(x)] <- ""
    vapply(x, function(s) {
        s <- gsub("[\u2018\u2019]", "'", s)
        s <- gsub("[\u201c\u201d]", "\"", s)
        lines <- strsplit(s, "\n", fixed = TRUE)[[1L]]
        ## timings, e.g. [35s/39s], [19m/19m], [1.2s]
        lines <- gsub("\\s*\\[[0-9.]+[smh](/[0-9.]+[smh])?\\]", "", lines)
        ## timing tables in examples/tests output
        drop <- grepl("with CPU \\(user \\+ system\\) or elapsed time", lines) |
            grepl("^\\s*user\\s+system\\s+elapsed\\s*$", lines) |
            grepl("^\\s*\\S+\\s+[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s*$", lines)
        lines <- lines[!drop]
        lines <- gsub("/tmp/Rtmp[[:alnum:]]+", "<TMPDIR>", lines)
        lines <- gsub("\\bfile[0-9a-f]{6,}\\b", "<FILE>", lines)
        lines <- gsub("\\S*/([^/[:space:]]+\\.Rcheck)/", "\\1/", lines)
        lines <- sub("\\s+$", "", lines)
        paste(lines[nzchar(lines)], collapse = "\n")
    }, "", USE.NAMES = FALSE)
}

## lines of b not in a, keeping b's order (multiset difference, so a line
## repeated more often in b than in a is reported)
line_diff <- function(a, b) {
    a <- if (nzchar(a)) strsplit(a, "\n", fixed = TRUE)[[1L]] else character()
    b <- if (nzchar(b)) strsplit(b, "\n", fixed = TRUE)[[1L]] else character()
    keep <- logical(length(b))
    for (i in seq_along(b)) {
        j <- match(b[i], a)
        if (is.na(j)) keep[i] <- TRUE else a <- a[-j]
    }
    b[keep]
}

compare_details <- function(old, new, include_better = FALSE) {
    packages <- intersect(old$Package, new$Package)
    db <- merge(old[old$Package %in% packages, ],
                new[new$Package %in% packages, ],
                by = c("Package", "Check"), all = TRUE,
                suffixes = c(".old", ".new"))
    ro <- status_rank(db$Status.old)
    rn <- status_rank(db$Status.new)
    ## only checks with an issue under at least one version
    db <- db[ro > 2L | rn > 2L, , drop = FALSE]
    ro <- status_rank(db$Status.old)
    rn <- status_rank(db$Status.new)
    no <- normalize_output(db$Output.old)
    nn <- normalize_output(db$Output.new)
    type <- ifelse(rn > ro, "worse",
            ifelse(rn < ro, "better",
            ifelse(no != nn, "changed", "same")))
    keep <- type %in% c("worse", "changed", if (include_better) "better")
    db <- db[keep, , drop = FALSE]
    no <- no[keep]; nn <- nn[keep]; type <- type[keep]
    st <- function(s) ifelse(is.na(s), "", s)
    res <- data.frame(Package   = db$Package,
                      Check     = db$Check,
                      Type      = type,
                      Old       = st(db$Status.old),
                      New       = st(db$Status.new),
                      stringsAsFactors = FALSE)
    res$Removed   <- Map(line_diff, nn, no, USE.NAMES = FALSE)
    res$Added     <- Map(line_diff, no, nn, USE.NAMES = FALSE)
    res$OldOutput <- no
    res$NewOutput <- nn
    ## flag packages whose own version differs between the two runs
    vo <- db$Version.old; vn <- db$Version.new
    ind <- !is.na(vo) & !is.na(vn) & vo != vn
    res$Package[ind] <- sprintf("%s [old version: %s, new version: %s]",
                                res$Package[ind], vo[ind], vn[ind])
    ord <- order(match(res$Type, c("worse", "changed", "better")),
                 res$Package, res$Check)
    res <- res[ord, , drop = FALSE]
    rownames(res) <- NULL
    attr(res, "only_old") <- setdiff(old$Package, new$Package)
    attr(res, "only_new") <- setdiff(new$Package, old$Package)
    res
}

format_changes <- function(res, full = FALSE, max_lines = 30L) {
    indent <- function(x, prefix) {
        if (!length(x)) return(character())
        n <- length(x)
        if (max_lines > 0L && n > max_lines)
            c(sprintf("%s[... %d earlier lines omitted]", prefix, n - max_lines),
              paste0(prefix, x[(n - max_lines + 1L):n]))
        else paste0(prefix, x)
    }
    fmt_out <- function(s) if (nzchar(s)) paste0("    ", strsplit(s, "\n", fixed = TRUE)[[1L]]) else "    (none)"
    hdr <- c(sprintf("Checks with changes: %d worse, %d changed messages, %d better",
                     sum(res$Type == "worse"), sum(res$Type == "changed"),
                     sum(res$Type == "better")),
             sprintf("Packages: %s",
                     paste(unique(sub(" .*", "", res$Package)), collapse = ", ")))
    for (w in c("only_old", "only_new"))
        if (length(p <- attr(res, w)))
            hdr <- c(hdr, sprintf("Results only in %s: %s",
                                  sub("only_", "", w), paste(p, collapse = ", ")))
    body <- lapply(seq_len(nrow(res)), function(i) {
        r <- res[i, ]
        c("",
          sprintf("Package: %s", r$Package),
          sprintf("Check: %s", r$Check),
          sprintf("Change: %s (%s -> %s)", r$Type,
                  if (nzchar(r$Old)) r$Old else "absent",
                  if (nzchar(r$New)) r$New else "absent"),
          if (full)
              c("  Old output:", fmt_out(r$OldOutput),
                "  New output:", fmt_out(r$NewOutput))
          else
              c(indent(r$Removed[[1L]], "  - "),
                indent(r$Added[[1L]],   "  + ")))
    })
    c(hdr, unlist(body))
}

checkChangesDetailed <-
function (args) {
    stopifnot(is.character(args))
    args.prefix <- sub("^(--.*?=)(.*)$", "\\1", args)
    args.suffix <- sub("^(--.*?=)(.*)$", "\\2", args)

    preclean <- TRUE
    include_better <- FALSE
    full <- FALSE
    max_lines <- 30L
    olddir <- newdir <- outdir <- NULL

    for (i in seq_along(args))
    switch (args.prefix[[i]],
            "--preclean" =
                preclean <- TRUE,
            "--no-preclean" =
                preclean <- FALSE,
            "--include-better" =
                include_better <- TRUE,
            "--full" =
                full <- TRUE,
            "--max-lines=" =
                max_lines <- as.integer(args.suffix[[i]]),
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
        outdir <- sprintf("CHANGES_DETAILED_old=%s_new=%s",
                          basename(olddir), basename(newdir))
    if (!dir.exists(outdir))
        dir.create(outdir)
    else if (preclean)
        unlink(file.path(outdir, c("CHANGES_detailed",
                                   "CHANGES_detailed.rds",
                                   "*_00install.out",
                                   "*_00check.log")))

    old <- tools::check_packages_in_dir_details(olddir, drop_ok = FALSE)
    new <- tools::check_packages_in_dir_details(newdir, drop_ok = FALSE)
    res <- compare_details(old, new, include_better = include_better)

    writeLines(out <- format_changes(res, full = full, max_lines = max_lines),
               file.path(outdir, "CHANGES_detailed"))
    writeLines(out)
    saveRDS(res, file = file.path(outdir, "CHANGES_detailed.rds"))

    package <- unique(sub(" .*", "", res$Package[res$Type != "better"]))
    for (zz in c("00install.out", "00check.log"))
    file.copy(file.path(newdir,
                        sprintf("rdepends_%s.Rcheck", package),
                        zz),
              file.path(outdir,
                        sprintf("%s_%s", package, zz)),
              overwrite = TRUE)

    invisible(res)
}

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    ch <- checkChangesDetailed(args)
}
