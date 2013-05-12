library("testthat")
library("lme4")

context("anova")
test_that("lmer", {
    fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
    fm0 <- update(fm1,.~.-Days)
    expect_that(anova(fm0,fm1),                        is_a("anova"))
    expect_warning(do.call(anova,list(fm0,fm1)),"assigning generic names")
})

