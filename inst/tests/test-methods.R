library("testthat")
library("lme4")

context("anova")
test_that("lmer", {
    fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
    fm0 <- update(fm1,.~.-Days)
    expect_that(anova(fm0,fm1),                        is_a("anova"))
    expect_warning(do.call(anova,list(fm0,fm1)),"assigning generic names")
})

context("bootMer")
test_that("bootMer", {
    ## testing bug-fix for ordering of sd/cor components in sd/cor matrix with >2 rows
    m1 <- lmer(strength~1+(cask|batch),Pastes)
    bb <- suppressWarnings(confint(m1,method="boot",nsim=3,quiet=TRUE))
    corvals <- bb[grep("^cor_",rownames(bb)),]
    expect_true(all(abs(corvals)<=1))
})
