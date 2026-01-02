#context("testing '||' notation for independent ranefs")

form1 <- Reaction ~ Days + (Days || Subject)
form2 <- Reaction ~ Days + (1 | Subject) + (0 + Days | Subject)
m1 <- lmer(form1, sleepstudy)
m2 <- lmer(form2, sleepstudy)

test_that("basic intercept + slope '||' works", {
  if (!getOption("lme4.doublevert.default", "split") == "split") skip()
  
  expect_equivalent(lFormula(form1, sleepstudy)$reTrms,
                    lFormula(form2, sleepstudy)$reTrms)
  expect_equivalent(fitted(m1),
                    fitted(m2))


  if (FALSE) {
    ## library(glmmTMB)
    m3 <- glmmTMB(form1, sleepstudy, REML = TRUE,
                  start = list(theta = log(c(25, 6))))
    dput(head(fitted(m3), 4L))
    }
    expect_equivalent(head(fitted(m2), 4L),
                      ## fitted values from glmmTMB:
                      c(252.917801422647, 272.708576602026,
                        292.499351781404, 312.290126960783),
                      tolerance = 5e-5)
})
