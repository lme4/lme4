## A model fitted inside a function stores 'control = ctrl' in its call,
## where 'ctrl' is local to that function and gone once it returns.
## fmF's formula is written in the fitting function, so the formula's
## environment holds 'ctrl'; fmG's formula is written here, where it does not.
fitDyestuff <- function() {
    ctrl <- lmerControl(optimizer = "bobyqa")
    lmer(Yield ~ 1|Batch, Dyestuff, control = ctrl)
}
fitWith <- function(form) {
    ctrl <- lmerControl(optimizer = "bobyqa")
    lmer(form, Dyestuff, control = ctrl)
}
fmF <- fitDyestuff()
fmG <- fitWith(Yield ~ 1|Batch)
ctrl0 <- lmerControl(optimizer = "bobyqa")
fm0 <- lmer(Yield ~ 1|Batch, Dyestuff, control = ctrl0)
## a gradient tolerance no fit can meet, with failure an error, makes
## every bootstrap refit fail: evidence that this control reached refit()
ctrlStrict <- lmerControl(optimizer = "bobyqa",
                          check.conv.grad = .makeCC("stop", tol = 1e-20))

test_that("bootMer() finds a local control object in the formula's environment", {
    b0 <- bootMer(fm0, fixef, nsim = 3, seed = 101)
    bF <- bootMer(fmF, fixef, nsim = 3, seed = 101)
    expect_equal(attr(bF, "bootFail"), 0)
    expect_equal(bF$t, b0$t)
})

test_that("bootMer() reports the caller-frame error when the formula's environment fails too", {
    expect_error(bootMer(fmG, fixef, nsim = 2), "'ctrl' not found")
})

test_that("bootMer() evaluates the stored control in 'call.env'", {
    b0 <- bootMer(fm0, fixef, nsim = 3, seed = 101)
    bList <- bootMer(fmG, fixef, nsim = 3, seed = 101,
                     call.env = list(ctrl = ctrl0))
    bEnv <- bootMer(fmG, fixef, nsim = 3, seed = 101,
                    call.env = list2env(list(ctrl = ctrl0)))
    expect_equal(attr(bList, "bootFail"), 0)
    expect_equal(bList$t, b0$t)
    expect_equal(bEnv$t, b0$t)
})

test_that("a list 'call.env' falls back to the caller's frame", {
    ctrl <- ctrl0
    bList <- bootMer(fmG, fixef, nsim = 3, seed = 101, call.env = list())
    b0 <- bootMer(fm0, fixef, nsim = 3, seed = 101)
    expect_equal(bList$t, b0$t)
})

test_that("'call.env' takes precedence over the formula's environment", {
    bStrict <- suppressWarnings(
        bootMer(fmF, fixef, nsim = 2, seed = 101,
                call.env = list(ctrl = ctrlStrict)))
    expect_equal(attr(bStrict, "bootFail"), 2)
    expect_true(all(is.na(bStrict$t)))
})

test_that("by default the caller's frame takes precedence over the formula's environment", {
    ctrl <- ctrlStrict
    bStrict <- suppressWarnings(bootMer(fmF, fixef, nsim = 2, seed = 101))
    expect_equal(attr(bStrict, "bootFail"), 2)
})

test_that("confint(method = \"boot\") passes 'call.env' to bootMer()", {
    ci <- suppressWarnings(
        confint(fmG, method = "boot", nsim = 3, seed = 101, quiet = TRUE,
                call.env = list(ctrl = ctrl0)))
    expect_equal(dim(ci), c(3L, 2L))
    expect_false(anyNA(ci))
})

test_that("bootMer() rejects a 'call.env' that is neither a list nor an environment", {
    expect_error(bootMer(fm0, fixef, nsim = 2, call.env = "ctrl"),
                 "must be an environment or a list")
})
