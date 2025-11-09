## transient patches for moving formula manipulation machinery to lme4 without breaking downstream packages

mkWarnFun <- function(FUN) {
  fn <- function(...) {
    msg <- sprintf("the %s function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", sQuote(FUN))
    rlang::warn(msg, .frequency = "once", .frequency_id = FUN)
    reformulas_fun <- getExportedValue("reformulas", FUN)
    reformulas_fun(...)
  }
  assign(FUN, fn, envir = parent.frame())
}

for (f in c("findbars","subbars", "nobars",
            "mkReTrms", "expandDoubleVerts", "isNested")) {
  mkWarnFun(f)
}

