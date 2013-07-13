library(lme4)
library(testthat)
d <- data.frame(x=runif(100),f=factor(rep(1:10,10)))
set.seed(101)
u <- rnorm(10)
d <- transform(d,y=rnorm(100,1+2*x+u[f],0.2))
d0 <- d
d[c(3,5,7),"x"] <- NA

## 'omit' and 'exclude' are the only choices under which
##  we will see NA values in the results
fm0 <- lmer(y~x+(1|f),data=d0)
## no 'na.action' attribute because no NAs in this data set
stopifnot(is.null(attr(model.frame(fm0),"na.action")))
fm1 <- update(fm0,data=d)
## no NAs in predict or residuals because na.omit
stopifnot(!any(is.na(predict(fm1))))
stopifnot(!any(is.na(residuals(fm1))))
fm2 <- update(fm1,na.action="na.exclude")
## no NAs in predict or residuals because na.omit
nNA <- sum(is.na(d$x))
stopifnot(sum(is.na(predict(fm2)))==nNA)
stopifnot(sum(is.na(residuals(fm2)))==nNA)
expect_error(fm3 <- lmer(y~x+(1|f),data=d,na.action="na.pass"),
             "infinite or missing values")

refit(fm0)
refit(fm1)
refit(fm2)

refit(fm0,runif(100))
refit(fm1,runif(100))
refit(fm2,runif(100))
