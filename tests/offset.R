## simple examples with offsets, to exercise methods etc.

library(lme4)
## generate a basic Gamma/random effects sim
set.seed(101)
d <- expand.grid(block=LETTERS[1:26],rep=1:100)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4+3*d$x  ## version without random effects
d$eta <- d$eta0+reff_f[d$block]

## lmer() test:
d$mu <- d$eta
d$y <- rnorm(nrow(d),mean=d$mu,sd=1)

fm1 <- lmer(y~x+(1|block),data=d)
fm1off <- lmer(y~x+(1|block)+offset(3*x),data=d)

## check equality
stopifnot(all.equal(fixef(fm1)[2]-3,fixef(fm1off)[2]))

p0 <- predict(fm1)
p1 <- predict(fm1,newdata=d)
p2 <- predict(fm1off,newdata=d)
stopifnot(all.equal(p0,p1,p2))


## glmer() test:
d$mu <- exp(d$eta)
d$y <- rpois(nrow(d),d$mu)

gm1 <- glmer(y~x+(1|block),data=d,family=poisson,
             control=glmerControl(check.conv.grad="ignore"))
gm1off <- glmer(y~x+(1|block)+offset(3*x),data=d,family=poisson,
                control=glmerControl(check.conv.grad="ignore"))

## check equality
stopifnot(all.equal(fixef(gm1)[2]-3,fixef(gm1off)[2],tolerance=3e-4))

p0 <- predict(gm1)
p1 <- predict(gm1,newdata=d)
p2 <- predict(gm1off,newdata=d)
stopifnot(all.equal(p0,p1,p2))

## FIXME: should also test simulations
