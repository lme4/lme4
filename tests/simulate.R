library(lme4)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
s1 <- simulate(fm1,seed=101)[[1]]
s2 <- simulate(fm1,seed=101,use.u=TRUE)

## binomial (2-column and prob/weights)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
gm2 <- glmer(incidence/size ~ period + (1 | herd), weights=size,
              data = cbpp, family = binomial)

s1 <- simulate(gm1,seed=101)[[1]]
s2 <- simulate(gm2,seed=101)[[1]]
stopifnot(all.equal(s1[,1]/rowSums(s1),s2))
s3 <- simulate(gm1,seed=101,use.u=TRUE)

## binomial (factor): Kubovy bug report 1 Aug 2013
d <- data.frame(y=factor(rep(letters[1:2],each=100)),
                f=factor(rep(1:10,10)))
g1 <- glmer(y~(1|f),data=d,family=binomial)
simulate(g1,nsim=10)

## test explicitly stated link function
gm3 <- glmer(cbind(incidence, size - incidence) ~ period +
             (1 | herd), data = cbpp, family = binomial(link="logit"))
s4 <- simulate(gm3,seed=101)[[1]]
stopifnot(all.equal(s1,s4))

cbpp$obs <- factor(seq(nrow(cbpp)))
gm4 <- glmer(cbind(incidence, size - incidence) ~ period +
             (1 | herd) + (1|obs), data = cbpp, family = binomial)

s5 <- simulate(gm4,seed=101)[[1]]
s6 <- simulate(gm4,seed=101,use.u=TRUE)[[1]]

## Bernoulli
## works, but too slow
if (FALSE) {
  data(guImmun,package="mlmRev")
  g1 <- glmer(immun~kid2p+mom25p+ord+ethn+momEd+husEd+momWork+rural+pcInd81+
              (1|comm/mom),family="binomial",data=guImmun)
  s2 <- simulate(g1)
}

set.seed(101)
d <- data.frame(f=rep(LETTERS[1:10],each=10))
d$x <- runif(nrow(d))
u <- rnorm(10)
d$eta <- with(d,1+2*x+u[f])
d$y <- rbinom(nrow(d),plogis(d$eta),size=1)

g1 <- glmer(y~x+(1|f),data=d,family="binomial")
## tolPwrss=1e-5: no longer necessary

if (FALSE) {
  allcoef <- function(x) {
    c(deviance(x),getME(x,"theta"),getME(x,"beta"))
  }
  tfun <- function(t) {
    gg <- try(glmer(y~x+(1|f),data=d,family="binomial",
                    control=glmerControl(tolPwrss=10^t)))
    if (inherits(gg,"try-error")) rep(NA,4) else allcoef(gg)
  }
  tvec <- seq(-4,-16,by=-0.25)
  tres <- cbind(tvec,t(sapply(tvec,tfun)))
}

s1 <- simulate(g1,seed=102)[[1]]

d$y <- factor(c("N","Y")[d$y+1])
g1B <- glmer(y~x+(1|f),data=d,family="binomial") ## ,tolPwrss=1e-5)
s1B <- simulate(g1B,seed=102)[[1]]
stopifnot(all.equal(s1,as.numeric(s1B)-1))

## another Bernoulli
data(Contraception,package="mlmRev")
fm1 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
s3 <- simulate(fm1)

d$y <- rpois(nrow(d),exp(d$eta))
g2 <- glmer(y~x+(1|f),data=d,family="poisson")
s4 <- simulate(g2)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
s5 <- simulate(fm1)
