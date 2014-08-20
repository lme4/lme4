## form <- y ~ x + (x | g) + d(~(x | g))
## formSplit <- lme4:::splitregen(form, "d")[1:2]
## stopifnot(do.call(identical, unname(lapply(formSplit, environment))))
