## Tests to prevent accidentally committing debugging statements.

lme4_find_src_dir <- function() {
  pkg_dir <- find.package("lme4", quiet = TRUE)
  src_dir <- file.path(pkg_dir, "src")
  ## C++ source files are only present in the source tree, not the installed package
  if (length(list.files(src_dir, pattern = "\\.cpp$")) > 0L) src_dir else NULL
}

lme4_find_r_dir <- function() {
  pkg_dir <- find.package("lme4", quiet = TRUE)
  r_dir <- file.path(pkg_dir, "R")
  if (length(list.files(r_dir, pattern = "\\.R$")) > 0L) r_dir else NULL
}

test_that("no active debug flags in C++ source files", {
  src_dir <- lme4_find_src_dir()
  skip_if(is.null(src_dir),
          "C++ source files not available (skipping in installed package)")
  
  cpp_files <- list.files(src_dir, pattern = "\\.(cpp|h|hpp)$",
                          full.names = TRUE)
  expect_gt(length(cpp_files), 0L)   ## sanity check: we found files
  
  bad_lines <- character(0L)
  for (f in cpp_files) {
    lines <- readLines(f, warn = FALSE)
    ## Flag any line where a 'debug' variable is set to a non-zero value.
    ## Catches:  int debug=1;   debug = 2;   debug=10;   etc.
    ## Ignores:  int debug=0;   (safe — debug disabled)
    ##           // debug=1;    (single-line comment)
    ##           /* debug=1; */ (block comment — line starts with '/', so
    ##                           the leading [^/]* anchor prevents a match)
    idx <- grep("^[^/]*\\bdebug\\s*=\\s*[1-9][0-9]*", lines)
    if (length(idx) > 0L) {
      bad_lines <- c(bad_lines,
                     sprintf("%s:%d: %s",
                             basename(f), idx, trimws(lines[idx])))
    }
  }
  
  expect_length(bad_lines, 0L)
})

test_that("lmer and glmer do not write to stdout", {
    out <- capture.output(
        suppressMessages(
            suppressWarnings({
                fit_lmer <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
                fit_glmer <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                                   data = cbpp, family = binomial)
                invisible(list(fit_lmer, fit_glmer))
            })
        )
    )
    expect_length(out, 0L)
})

test_that("no active debugging calls in R source files", {
    r_dir <- lme4_find_r_dir()
    skip_if(is.null(r_dir),
            "R source files not available (skipping in installed package)")

    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    expect_gt(length(r_files), 0L)   ## sanity check: we found files

    bad_lines <- character(0L)
    for (f in r_files) {
        lines <- readLines(f, warn = FALSE)
        ## Flag browser() calls that are not commented out.
        ## Catches:  browser()   browser(expr = ...)
        idx <- grep("\\bbrowser\\s*\\(", lines)
        if (length(idx) > 0L) {
            bad_lines <- c(bad_lines,
                           sprintf("%s:%d: %s",
                                   basename(f), idx, trimws(lines[idx])))
        }
    }

    expect_length(bad_lines, 0L)
})
