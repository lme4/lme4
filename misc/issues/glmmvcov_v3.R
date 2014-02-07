## my (RHC) installation of gitHub lme4:
##library(lme4, lib.loc="/Users/rhbc/Documents/Misc/tmp_Rlibs/")
library(lme4)
## install.packages("ordinal", repos="http://R-Forge.R-project.org")
stopifnot(require(ordinal))
## Get example data:
data(wine, package="ordinal")
wine <- within(wine,
               binresp <- factor(response < 50, labels=c("high", "low"))
               )
## head(wine)

## Function to compute hessian and gradient:
deriv12 <- function(fun, x, delta=1e-4, fx=NULL, ...) {
### Compute gradient and Hessian at the same time (to save computing
### time)
    nx <- length(x)
    fx <- if(!is.null(fx)) fx else fun(x, ...)
    stopifnot(is.numeric(fx))
    stopifnot(length(fx) == 1)
    H <- array(NA, dim=c(nx, nx))
    g <- numeric(nx)
    for(j in 1:nx) {
        ## Diagonal elements:
        xadd <- xsub <- x
        xadd[j] <- x[j] + delta
        xsub[j] <- x[j] - delta
        fadd <- fun(xadd, ...)
        fsub <- fun(xsub, ...)
        H[j, j] <- (fadd - 2 * fx + fsub) / delta^2
        g[j] <- (fadd - fsub) / (2 * delta)
        ## Off diagonal elements:
        for(i in 1:nx) {
            if(i >= j) break
            xaa <- xas <- xsa <- xss <- x
            xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
            xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
            xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
            xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
            H[i, j] <- H[j, i] <-
                (fun(xaa, ...) - fun(xas, ...) -
                 fun(xsa, ...) + fun(xss, ...)) /
                     (4 * delta^2)
        }
    }
    list(gradient = g, Hessian = H)
}


## Variance-covariance matrix / standard errors from glmer do not
## (always) correspond to the curvature in the log-likelihood
## function:

## wine data example:
gm2 <- glmer(binresp ~ temp + contact + (1|judge), data=wine,
             family = binomial)
## Compute hessian for glmer fit:
gm2Devfun <- update(gm2,devFunOnly=TRUE)
parVals <- c(gm2@theta, fixef(gm2))
d12 <- deriv12(fun=gm2Devfun, x=parVals)
## gradient:
d12$gradient
## Hessian (of the log-lik wrt. all parameters):
gm2.hess <- d12$Hessian/2 ## dev -> log-lik scale.
## Standard errors
(se.hess <- sqrt(diag(solve(gm2.hess)))[-1]) ## from Hessian
(se.glmer <- sqrt(diag(vcov(gm2)))) ## from glmer
se.glmer - se.hess
abs(se.glmer - se.hess)/se.hess * 100
## around 10% relative difference.
## Try with ordinal::clmm fit:
fm2 <- clmm(binresp ~ temp + contact + (1|judge), data=wine)
(se.clmm <- sqrt(diag(vcov(fm2)))[-4])
abs(se.clmm - se.hess)/se.hess * 100
## < .01% relative difference.

## cbpp example:
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
## Compute hessian:
gm1Devfun <- update(gm1,devFunOnly=TRUE)
parVals <- c(gm1@theta, fixef(gm1))
d12 <- deriv12(fun=gm1Devfun, x=parVals)
## gradient:
d12$gradient
## Hessian (of the log-lik wrt. all parameters):
gm1.hess <- d12$Hessian/2 ## dev -> log-lik scale.
## Standard errors
(se.hess <- sqrt(diag(solve(gm1.hess)))[-1]) ## from Hessian
(se.glmer <- sqrt(diag(vcov(gm1)))) ## from glmer
se.glmer - se.hess
abs(se.glmer - se.hess)/se.hess * 100
## around 1% relative difference.

## I think the small difference for the cbpp fit compared to the wine
## fit is due to 1) much more information and
## 2) less correlation between theta and the other parameters
## (probably affected by 1)).

## ad.1): There are more than 10 times as many Bernoulli trials in the
## cbpp data compared to the wine data:
with(cbpp, sum(c(incidence, size)))
nrow(wine)

## ad.2) The correlation with theta is much smaller in the cbpp fit:
round(cov2cor(solve(gm1.hess)), 3) ## cbpp
round(cov2cor(solve(gm2.hess)), 3) ## wine

sessionInfo()

