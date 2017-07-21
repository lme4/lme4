stopifnot(require("testthat"), require("lme4"))

set.seed(101)
n <- 500
d <- data.frame(x=rnorm(n),
                f=factor(sample(1:10,n,replace=TRUE),
                         labels=LETTERS[1:10]),
                g=factor(sample(1:25,n,replace=TRUE),
                         labels=letters[1:25]))
d$y <- suppressMessages(simulate(~1+x+(1|f)+(x|g),family=binomial,
                                 newdata=d,
                                 newparams=list(beta=c(0,1),
                                                theta=c(1,1,2,1)))[[1]])
fm1 <- glmer(y~(1|f)+(x|g),family=binomial,data=d)

test_that("ranef", {
    expect_warning(ranef(fm1,transf=exp),"additional arguments")
})
test_that("dotplot_ranef", {
    rr <- ranef(fm1,condVar=TRUE)
    expect_is(lattice::dotplot(rr,scales=list(x = list(relation = 'free')))$g,
              "trellis")
    expect_is(lattice::dotplot(rr,transf=exp,
                     scales=list(x = list(relation = 'free')))$g,
              "trellis")
    expect_is(as.data.frame(rr),"data.frame")
    rr0 <- ranef(fm1)
    expect_is(as.data.frame(rr0),"data.frame")
})
