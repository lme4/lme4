library("testthat")
library("lme4")

context("factor handling in grouping variables")

test_that("factors", {

    set.seed(101)
    d <- data.frame(x=runif(1000),y=runif(1000),f1=rep(1:10,each=100),f2=rep(1:10,100))
    d2 <- transform(d,f1=factor(f1),f2=factor(f2))
    expect_that(lm1 <- lmer(y~x+(1|f1/f2),data=d), is_a("lmerMod"))
    expect_that(lm2 <- lmer(y~x+(1|f1/f2),data=d2),is_a("lmerMod"))
    expect_equivalent(lm1,lm2)

})
