library(lme4)
devtools::load_all() ## get version with na.action.merMod)
library(testthat)
d <- lme4::sleepstudy
d$Reaction[1] <- NA
fm1 <- lmer(Reaction ~ Days + (Days | Subject), d)
## correct dimensions with na.exclude as well ?
fm2 <- update(fm1, na.action = na.exclude)

expect_equal(length(residuals(fm1)), nrow(model.frame(fm1)))
expect_equal(length(residuals(fm2)), nrow(d))
na.action(fm2)
w <- weights(fm2)  ## simple numeric vector
length(weights(fm2))

weighted.residuals <- function (obj, drop0 = TRUE)  {
    na.act <- na.action(obj)
    w <- naresid(na.act, weights(obj))
    ## current version doesn't work if residuals() already uses napredict() ?
    ## (double-correction)
    ## what do we do to get working residuals??
    r <- residuals(obj) ## naresid(na.act, residuals(obj))
    if (!is.null(w)) 
        r <- r * sqrt(w)
    if (inherits(obj, "glm")) 
        w <- weights(obj, "prior")
    if (drop0 && !is.null(w)) {
        if (is.matrix(r)) 
            r[w != 0, , drop = FALSE]
        else r[w != 0]
    }
    else r
}

length(weights(fm2))
w1 <- try(stats::weighted.residuals(fm2))
w2 <- weighted.residuals(fm2)



## TO DO:
## tests in lme4
## test some base-R objects with na.exclude()

d2 <- mtcars
d2$mpg[1] <- NA
m0 <- lm(mpg ~ hp, data = d2)
stats::weighted.residuals(m0)

length(weighted.residuals(m0))
na.act <- na.action(m0)
m1 <- update(m0, na.action = na.exclude)
stopifnot(length(weighted.residuals(m0)) == nrow(na.omit(d2)))
stopifnot(length(weighted.residuals(m1)) == nrow(d2))
