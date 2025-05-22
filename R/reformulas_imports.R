## transient patches for moving formula manipultion machinery to lme4 without breaking downstream packages
subbars <- function(...) {
  rlang::warn("the subbars() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "subbars")
  eval.parent(reformulas::subbars(...))
}

nobars <- function(...) {
  rlang::warn("the nobars() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "nobars")
  eval.parent(reformulas::nobars(...))
}

findbars <- function(...) {
  rlang::warn("the findbars() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "findbars")
  eval.parent(reformulas::findbars(...))
}

mkReTrms <- function(...) {
  rlang::warn("the mkReTrms() function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainter to do so.", .frequency = "once", .frequency_id = "mkReTrms")
  eval.parent(reformulas::mkReTrms(...))
}

