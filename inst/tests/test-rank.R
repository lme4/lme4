library("testthat")
library("lme4")

context("testing fixed-effect design matrices for full rank")

test_that("lmerRank", {
    set.seed(101)
    n <- 20
    x <- y <- rnorm(n)
    z <- rnorm(n)
    r <- sample(1:5, size=n, replace=TRUE)
    d <- data.frame(x,y,z,r)
    d$y2 <- d$y + c(0.001,rep(0,n-1))

    expect_error(lmer( z ~ x + y + (1|r), data=d),"rank of X")
    ## should work:
    expect_that(lmer( z ~ x + y2 + (1|r), data=d), is_a("lmerMod"))

    d2 <- expand.grid(a=factor(1:4),b=factor(1:4),rep=1:10)
    n <- nrow(d2)
    d2 <- transform(d2,r=sample(1:5, size=n, replace=TRUE),
                     z=rnorm(n))
    d2 <- subset(d2,!(a=="4" & b=="4"))
    expect_error(lmer( z ~ a*b + (1|r), data=d2),"rank of X")
    d2 <- transform(d2, ab=droplevels(interaction(a,b)))
    ## should work:
    expect_that(lmer( z ~ ab + (1|r), data=d2), is_a("lmerMod"))
})

test_that("glmerRank", {
    set.seed(101)
    n <- 100
    x <- y <- rnorm(n)
    z <- rbinom(n,size=1,prob=0.5)
    r <- sample(1:5, size=n, replace=TRUE)
    d <- data.frame(x,y,z,r)
    ## d$y2 <- d$y + c(0.001,rep(0,n-1))  ## too small: get convergence failures
    ## FIXME: figure out how small a difference will still fail?
    d$y2 <- rnorm(n)

    expect_error(glmer( z ~ x + y + (1|r), data=d, family=binomial), "rank of X")
    expect_that(glmer( z ~ x + y2 + (1|r), data=d, family=binomial),is_a("glmerMod"))
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



