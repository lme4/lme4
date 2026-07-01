## transient patches for moving formula manipulation machinery to lme4 without breaking downstream packages

## tracks which one-time deprecation warnings have already fired this session
.reformulas_warned <- new.env(parent = emptyenv())

warn_once <- function(msg, id) {
  if (!exists(id, envir = .reformulas_warned, inherits = FALSE)) {
    assign(id, TRUE, envir = .reformulas_warned)
    warning(msg, call. = FALSE)
  }
}

mkWarnFun <- function(FUN) {
  fn <- function(...) {
    msg <- sprintf("the %s function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainer to do so.", sQuote(FUN))
    warn_once(msg, FUN)
    reformulas_fun <- getExportedValue("reformulas", FUN)
    reformulas_fun(...)
  }
  assign(FUN, fn, envir = parent.frame())
}

for (f in c("findbars","subbars", "nobars",
            "mkReTrms", "expandDoubleVerts", "isNested")) {
  mkWarnFun(f)
}

## handle reOnly separately, it needs to be wrapped to return a formula properly
## reformulas 0.4.5 will fix this (can then move back to the deprecation block above)
reOnly <- function(...) {
  res <- lme4_reOnly(...)
}
