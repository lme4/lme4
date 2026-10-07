## A model fitted inside a function stores 'control = ctrl' in its call,
## where 'ctrl' is local to that function and gone once it returns.
fitDyestuff <- function() {
    ctrl <- lmerControl(optimizer = "bobyqa")
    lmer(Yield ~ 1|Batch, Dyestuff, control = ctrl)
}
fmF <- fitDyestuff()
ctrl0 <- lmerControl(optimizer = "bobyqa")
fm0 <- lmer(Yield ~ 1|Batch, Dyestuff, control = ctrl0)

test_that("bootMer() cannot find a control object local to the fitting function", {
    expect_error(bootMer(fmF, fixef, nsim = 2), "'ctrl' not found")
})

test_that("bootMer() evaluates the stored control in 'call.env'", {
    b0 <- bootMer(fm0, fixef, nsim = 3, seed = 101)
    bList <- bootMer(fmF, fixef, nsim = 3, seed = 101,
                     call.env = list(ctrl = ctrl0))
    bEnv <- bootMer(fmF, fixef, nsim = 3, seed = 101,
                    call.env = list2env(list(ctrl = ctrl0)))
    expect_equal(attr(bList, "bootFail"), 0)
    expect_equal(bList$t, b0$t)
    expect_equal(bEnv$t, b0$t)
})

test_that("a list 'call.env' falls back to the caller's frame", {
    ctrl <- ctrl0
    bList <- bootMer(fmF, fixef, nsim = 3, seed = 101, call.env = list())
    b0 <- bootMer(fm0, fixef, nsim = 3, seed = 101)
    expect_equal(bList$t, b0$t)
})

test_that("the control found in 'call.env' reaches refit()", {
    ## a gradient tolerance no fit can meet, with failure an error, makes
    ## every bootstrap refit fail
    ctrlStrict <- lmerControl(optimizer = "bobyqa",
                              check.conv.grad = .makeCC("stop", tol = 1e-20))
    bStrict <- suppressWarnings(
        bootMer(fmF, fixef, nsim = 2, seed = 101,
                call.env = list(ctrl = ctrlStrict)))
    expect_equal(attr(bStrict, "bootFail"), 2)
    expect_true(all(is.na(bStrict$t)))
})

test_that("confint(method = \"boot\") passes 'call.env' to bootMer()", {
    ci <- suppressWarnings(
        confint(fmF, method = "boot", nsim = 3, seed = 101, quiet = TRUE,
                call.env = list(ctrl = ctrl0)))
    expect_equal(dim(ci), c(3L, 2L))
    expect_false(anyNA(ci))
})

test_that("bootMer() rejects a 'call.env' that is neither a list nor an environment", {
    expect_error(bootMer(fm0, fixef, nsim = 2, call.env = "ctrl"),
                 "must be an environment or a list")
})
