library(lme4)

options(nwarnings = 1000)

if(getRversion() < "3.2.0") {
    if(interactive()) break # gives an error
    else q() # <- undesirable when interactive !
}

## Try all "standard" (statistical) S3 methods:
.S3generics <- function(class) {
    s3m <- .S3methods(class=class)
    ii <- attr(s3m, "info")
    ii[!ii[, "isS4"], "generic"]
}

set.seed(12)
d <- data.frame(
  g = sample(c("A","B","C","D","E"), 250, replace=TRUE),
  y1 = runif(250, max=100),
  y2 = sample(c(0,1), 250, replace=TRUE)
)

fm3.1 <- lmList(y1 ~ 1 | g, data=d)
fm3.2 <- lmList(y2 ~ 1 | g, data=d, family=binomial)

data(Orthodont, package="nlme")
Orthodont <- as.data.frame(Orthodont) # no "groupedData"
fm2 <- lmList(distance ~ age | Subject, Orthodont)

s3fn <- .S3generics(class= class(fm3.1)[1]) ## works for "old and new" class

noquote(s3fn <- s3fn[s3fn != "print"])# <-- it is show() not print() that works
## [1] coef    confint fitted fixef     formula logLik  pairs  plot
## [9] predict qqnorm  ranef  residuals sigma   summary update

## In lme4 1.1-7 (July 2014), only these worked:
##  coef(), confint(), formula(), logLik(), summary(), update()

## pairs() is excluded for fm3.1 which has only intercept:
## no errors otherwise:
evs <- sapply(s3fn[s3fn != "pairs"], do.call, args = list(fm3.1))

cls <- sapply(evs, function(.) class(.)[1])
clsOk <- cls[c("confint", "fixef", "formula", "logLik",
               "ranef", "sigma", "summary", "update")]
stopifnot(identical(unname(clsOk),
                    c("lmList4.confint", "numeric", "formula", "logLik",
                      "ranef.lmList", "numeric", "summary.lmList", "lmList4")))

## --- fm2 --- non-trivial X: can use pairs(), too:
evs2 <- sapply(s3fn, do.call, args = list(fm2))

## --- fm3.2 --- no failures for this "glmList" :
ss <- function(...) suppressMessages(suppressWarnings(...))
ss(evs3.2 <- sapply(s3fn[s3fn != "pairs"], do.call,
                                  args = list(fm3.2)))

## --- fm4 ---
evs4 <- sapply(s3fn, function(fn)
    tryCatch(do.call(fn, list(fm4)), error=function(e) e))
length(warnings())
summary(warnings()) ##  4 kinds;  glm.fit: fitted probabilities numerically 0 or 1 occurred

str(sapply(evs4, class)) # more errors than above
isok4 <- !sapply(evs4, is, class2="error")

## includes a nice pairs():
evs4[isok4]
## Error msgs of those with errors, first 5, now 3, then 2 :
str(errs4 <- lapply(evs4[!isok4], conditionMessage))
## $ logLik : chr "log-likelihood not available with NULL fits"
## $ summary: chr "subscript out of bounds"
stopifnot(length(errs4) <= 2)

