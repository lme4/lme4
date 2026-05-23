## see https://github.com/lme4/lme4/pull/811
## https://github.com/lme4/lme4/pull/954
##  (which raised the initial question of whether we could
##   move the assign() stuff from lFormula/glFormula into their
##   refactored/modular core, and whether we would need an additional
##   level of parent.frame)
## https://github.com/lme4/lme4/issues/481 (overwriting weights/offsets in global frame)

library(lme4)
library(testthat)

## example 1 from krlmr
my_weights <- rep(1, 6L)
formula <- mpg ~ disp + (1 | cyl)

lf1 <- local({
  my_weights = rep(2, 6L)
  lme4::lFormula(formula, data = head(mtcars), weights = my_weights)$fr[["(weights)"]]
})
print(lf1) ## matches global, not local (confusing, but defensible, as
           ##  it's looking in the environment of the formula)

rm(my_weights) ## clean up
rm(formula)

## example 2 from krlmr
formula = mpg ~ disp + (1 | cyl)

local({
  tbl = head(mtcars)
  try(lme4::lFormula(formula, data = tbl, weights = tbl$hp)$fr)
})

## fails (object 'tbl' not found)

## this works:
local({
  tbl = head(mtcars)
  lme4::lFormula(formula, data = tbl, weights = hp)$fr[["(weights)"]]
})


## with lFormula(...) this breaks model.frame already
## match.call is
##   lFormula(formula = ..1, data = ..2, weights = ..3)
## mf has turned into
##  stats::model.frame(data = ..2, weights = ..3, drop.unused.levels = TRUE, 
##         formula = mpg ~ disp + (1 + cyl))

get_comp <- function(..., w = c("(weights)", "(offset)")) {
    w <- match.arg(w)
    lf <- do.call(lFormula, list(...))
    ## lf <- lFormula(...)
    lf$fr[[w]]
    
}

lf1 <- local({
  my_weights = rep(2, 6L)
  get_comp(formula, data = head(mtcars), weights = my_weights)
})

lf1

## try some simulation tests
##  * new code doesn't break simulations
##  * we don't create new 'weights' vars in the formula environment,
##   or overwrite existing ones

cbpp$obs <- factor(seq(nrow(cbpp)))
fit_cbpp_0 <- glmer(cbind(incidence, size-incidence) ~ 1 + (1|herd),
                    cbpp, family=binomial)
## include fixed effect of period
fit_cbpp_1 <- update(fit_cbpp_0, . ~ . + period)
## include observation-level RE
fit_cbpp_2 <- update(fit_cbpp_1, . ~ . + (1|obs))

environment(formula(fit_cbpp_2))  ## global

chk_wts_absent <- function() {
    testthat::expect_identical(find("weights"), c(".GlobalEnv", "package:stats"))
}
weights <- wts_orig <- c(NA_real_, NA_real_)

expect_is(simulate(fit_cbpp_2), "data.frame")
p1 <- simulate(fit_cbpp_2, re.form = NULL, seed = 101)
sim <- simulate(formula(fit_cbpp_2), newdata = cbpp,
         family = binomial,
         newparams = list(theta = c(1,1),
                          beta = rep(0,4)))
chk_wts_absent()
expect_identical(weights, wts_orig)
