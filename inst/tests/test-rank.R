library("testthat")
library("lme4")

context("testing fixed-effect design matrices for full rank")

test_that("lmerRank", {
    set.seed(101)
    n <- 20
    x <- y <- rnorm(n)
    d <- data.frame(x,y,
		    z = rnorm(n),
		    r = sample(1:5, size=n, replace=TRUE),
		    y2 = y + c(0.001, rep(0,n-1)))
    expect_message(fm <- lmer( z ~ x + y + (1|r), data=d),
		   "fixed-effect model matrix is .*rank deficient")
    expect_equal(nrow(anova(fm)), 1L)
    expect_error(lmer( z ~ x + y + (1|r), data=d,
		      control=lmerControl(check.rankX="stop")),
		 "rank deficient")
    expect_error(lmer( z ~ x + y + (1|r), data=d,
		      control=lmerControl(check.rankX="ignore")),
		 "not positive definite")
    ## should work:
    expect_is(lmer( z ~ x + y2 + (1|r), data=d), "lmerMod")

    d2 <- expand.grid(a=factor(1:4),b=factor(1:4),rep=1:10)
    n <- nrow(d2)
    d2 <- transform(d2,r=sample(1:5, size=n, replace=TRUE),
                     z=rnorm(n))
    d2 <- subset(d2,!(a=="4" & b=="4"))
    expect_error(lmer( z ~ a*b + (1|r), data=d2,
                      control=lmerControl(check.rankX="stop")),
                 "rank deficient")
    expect_message(fm <- lmer( z ~ a*b + (1|r), data=d2),
                   "fixed-effect model matrix is rank deficient")
    d2 <- transform(d2, ab=droplevels(interaction(a,b)))
    ## should work:
    expect_is(fm2 <- lmer( z ~ ab + (1|r), data=d2), "lmerMod")
    expect_equal(logLik(fm), logLik(fm2))
    expect_equal(sum(anova(fm)[, "Df"]), anova(fm2)[, "Df"])
    expect_equal(sum(anova(fm)[, "Sum Sq"]), anova(fm2)[, "Sum Sq"])
})

test_that("glmerRank", {
    set.seed(111)
    n <- 100
    x <- y <- rnorm(n)
    d <- data.frame(x, y,
		    z = rbinom(n,size=1,prob=0.5),
		    r = sample(1:5, size=n, replace=TRUE),
                    y2 = ## y + c(0.001,rep(0,n-1)), ## too small: get convergence failures
                    ## FIXME: figure out how small a difference will still fail?
                    rnorm(n))
    expect_message(fm <- glmer( z ~ x + y + (1|r), data=d, family=binomial),
                   "fixed-effect model matrix is rank deficient")
    expect_error(glmer( z ~ x + y + (1|r), data=d, family=binomial,
                       control=glmerControl(check.rankX="stop")),
                 "rank deficient.*rank.X.")
    expect_is(glmer( z ~ x + y2 + (1|r), data=d, family=binomial), "glmerMod")
})

test_that("nlmerRank", {
    set.seed(101)
    n <- 1000
    nblock <- 15
    x <- abs(rnorm(n))
    y <- rnorm(n)
    z <- rnorm(n,mean=x^y)
    r <- sample(1:nblock, size=n, replace=TRUE)
    d <- data.frame(x,y,z,r)
    ## save("d","nlmerRank.RData")  ## see what's going on with difference in contexts

    fModel <- function(a,b)  (exp(a)*x)^(b*y)
    fModf <- deriv(body(fModel), namevec = c("a","b"),
                   func = fModel)

    fModel2 <- function(a,b,c)  (exp(a+c)*x)^(b*y)
    fModf2 <- deriv(body(fModel2), namevec = c("a","b","c"),
                   func = fModel2)

    ## should be OK: fails in test mode?
    nlmer(y ~ fModf(a,b) ~ a|r, d, start = c(a=1,b=1))

    ## FIXME: this doesn't get caught where I expected
    expect_error(nlmer(y ~ fModf2(a,b,c) ~ a|r, d, start = c(a=1,b=1,c=1)),"Downdated VtV")

})



