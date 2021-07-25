library("testthat")
library("lme4")


test_that("factors", {

    set.seed(101)
    d <- data.frame(x=runif(1000),y=runif(1000),f1=rep(1:10,each=100),f2=rep(1:10,100))
    d2 <- transform(d,f1=factor(f1),f2=factor(f2))
    expect_that(lm1 <- lmer(y~x+(1|f1/f2),data=d), is_a("lmerMod"))
    expect_that(lm2 <- lmer(y~x+(1|f1/f2),data=d2),is_a("lmerMod"))
    expect_equivalent(lm1,lm2)

})

## this will fail/take a long time unless we handle interactions carefully
test_that("savvy interactions", {
    dd <- data.frame(y = 1:10000, f1 = factor(1:10000), f2 = factor(1:10000))
    F1 <- lFormula(y ~ 1 + (1|f1/f2), data =dd,
             control = lmerControl(check.nobs.vs.nlev = "ignore",
                                   check.nobs.vs.nRE = "ignore"))
    expect_equal(dim(F1$reTrms$Zt), c(20000, 10000))
})

test_that("savvy factor level ordering", {

    check_f <- function(n = 200, frac = 0.7, fix_order = TRUE, check_order = TRUE) {
        dd <- expand.grid(f1 = seq(n), f2 = seq(n))
        dd <- within(dd, {
            f1 <- factor(f1, levels = sample(unique(f1)))
            f2 <- factor(f2, levels = sample(unique(f2)))
        })
        dd <- dd[sample(nrow(dd), size = round(frac*nrow(dd)), replace = FALSE), ]
        dd <- within(dd, {
            f12 <- f1:f2
            f12d <- droplevels(f12)
        })
        new_levels <- with(dd, levels(`%i%`(f1,f2, fix.order = fix_order)))
        ## don't want to pay the cost of checking if unneeded {for benchmarking}
        if (fix_order && check_order) { stopifnot(identical(levels(dd$f12d), new_levels)) }
        return(TRUE)
    }

    ## should fail within check_f() if levels don't match
    expect_true(check_f(), "'savvy' factor levels match brute-force version")
    
    ## library(microbenchmark)
    ## set.seed(101)
    ## m1 <- microbenchmark(check_f(fix_order = TRUE, check_order = FALSE),
    ## check_f(fix_order = FALSE))
})

