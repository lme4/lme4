library("lme4")
library("testthat")

set.seed(101)
dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=1:9, rep=1:15,
          KEEP.OUT.ATTRS=FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

## mimic glmer.nb protocol

test_that("ok with Poisson masking", {
    poisson <- NA
    ## use shortened version of data for speed ...
    m.base <- glmer.nb(y ~ f1 + (1|g), data=dd[1:200,])
    expect_is(m.base,"merMod")
})

context("testing glmer refit")

test_that("glmer refit", {
            ## basic Poisson fit
            m.base <- glmer(y ~ f1*f2 + (1|g), data=dd, family=poisson)
            expect_equal(m.base@beta,(m.base.r <- refit(m.base))@beta,
                         tolerance = 1e-5)

            th <- lme4:::est_theta(m.base,limit=20,eps=1e-4,trace=FALSE)
            th0 <- structure(0.482681268108477, SE = 0.0244825021248148)
            th1 <- structure(0.482681277470945)
            th2 <- 0.482681268108477
            th3 <- 0.4826813
            ## NB update with raw number
            m.numth1 <- update(m.base,family=negative.binomial(theta=0.4826813))
            expect_equal(m.numth1@beta,(m.numth1.r <- refit(m.numth1))@beta)

            ## strip NB value
            m.symth4 <- update(m.base,family=negative.binomial(theta=c(th)))
            expect_equal(m.symth4@beta,(m.symth4.r <- refit(m.symth4))@beta)

            ## IDENTICAL numeric value to case #1 above
            m.symth6 <- update(m.base,family=negative.binomial(theta=th3))
            expect_equal(m.symth6@beta,(m.symth6.r <- refit(m.symth6))@beta)

            ## standard NB update with computed theta from est_theta (incl SE attribute)
            m.symth <- update(m.base,family=negative.binomial(theta=th))
            expect_equal(m.symth@beta,(m.symth.r <- refit(m.symth))@beta)

            ## NB update with equivalent value
            m.symth2 <- update(m.base,family=negative.binomial(theta=th0))
            expect_equal(m.symth2@beta,(m.symth2.r <- refit(m.symth2))@beta) 

            ## NB update with theta value (stored as variable, no SE) only
            m.symth3 <- update(m.base,family=negative.binomial(theta=th1))
            expect_equal(m.symth3@beta,(m.symth3.r <- refit(m.symth3))@beta)

            ## strip NB value (off by 5e-16)
            m.symth5 <- update(m.base,family=negative.binomial(theta=th2))
            expect_equal(m.symth5@beta,(m.symth5.r <- refit(m.symth5))@beta)
        })

## GH #399
test_that("na_exclude", {
    dd1 <- dd[1:200,]
    dd1$f1[1:5] <- NA
    expect_error(glmer.nb(y ~ f1 + (1|g), data=dd1, na.action=na.fail),
                 "missing values in object")
    m1 <- glmer.nb(y ~ f1 + (1|g), data=dd1, na.action=na.omit)
    m2 <- glmer.nb(y ~ f1 + (1|g), data=dd1, na.action=na.exclude)
    expect_equal(fixef(m1),fixef(m1))
    expect_equal(length(predict(m2))-length(predict(m1)),5)
})

## GH 423
test_that("start_vals", {
    dd1 <- dd[1:200,]
    g1 <- glmer.nb(y ~ f1 + (1|g), data=dd1)
    g2 <- glmer.nb(y ~ f1 + (1|g), data=dd1,
                   initCtrl=list(theta=getME(g1,"glmer.nb.theta")))
    expect_equal(fixef(g1),fixef(g2),tol=1e-5)
})
