## transient patches for moving formula manipulation machinery to lme4 without breaking downstream packages
subbars <- function(...) {
  rlang::warn("the subbars() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "subbars")
  reformulas::subbars(...)
}

nobars <- function(...) {
  rlang::warn("the nobars() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "nobars")
  reformulas::nobars(...)
}

findbars <- function(...) {
  rlang::warn("the findbars() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "findbars")
  reformulas::findbars(...)
}

mkReTrms <- function(...) {
  rlang::warn("the mkReTrms() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "mkReTrms")
  reformulas::mkReTrms(...)
}

expandDoubleVerts <- function(...) {
  rlang::warn("the expandDoubleVerts() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "expandDoubleVerts")
  reformulas::expandDoubleVerts(...)
}

