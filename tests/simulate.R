library(lme4)
library(testthat)

(testLevel <- lme4:::testLevel())
L <- load(system.file("testdata/lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))
if (getRversion()>"3.0.0") {
    ## saved fits are not safe with old R versions

fm1 <- fit_sleepstudy_1

s1 <- simulate(fm1,seed=101)[[1]]
s2 <- simulate(fm1,seed=101,use.u=TRUE)
s3 <- simulate(fm1,seed=101,nsim=10)
s4 <- simulate(fm1,seed=101,use.u=TRUE,nsim=10)
stopifnot(length(s3)==10,all(sapply(s3,length)==180),
          length(s4)==10,all(sapply(s4,length)==180))

## binomial (2-column and prob/weights)
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_3

gm1_s1 <- simulate(gm1,seed=101)[[1]]
gm1_s2 <- simulate(gm2,seed=101)[[1]]
stopifnot(all.equal(gm1_s1[,1]/rowSums(gm1_s1),gm1_s2))
gm1_s3 <- simulate(gm1,seed=101,use.u=TRUE)
gm1_s4 <- simulate(gm1,seed=101,nsim=10)
gm1_s5 <- simulate(gm2,seed=101,nsim=10)
stopifnot(length(gm1_s4)==10,all(sapply(gm1_s4,ncol)==2),all(sapply(gm1_s4,nrow)==56))
stopifnot(length(gm1_s5)==10,all(sapply(gm1_s5,length)==56))

## binomial (factor): Kubovy bug report 1 Aug 2013
d <- data.frame(y=factor(rep(letters[1:2],each=100)),
                f=factor(rep(1:10,10)))
g1 <- glmer(y~(1|f),data=d,family=binomial)
s6 <- simulate(g1,nsim=10)
stopifnot(length(s6)==10,all(sapply(s6,length)==200))

## test explicitly stated link function
gm3 <- glmer(cbind(incidence, size - incidence) ~ period +
             (1 | herd), data = cbpp, family = binomial(link="logit"))
s4 <- simulate(gm3,seed=101)[[1]]
stopifnot(all.equal(gm1_s1,s4))

cbpp$obs <- factor(seq(nrow(cbpp)))
gm4 <- fit_cbpp_2
## glmer(cbind(incidence, size - incidence) ~ period +
##             (1 | herd) + (1|obs), data = cbpp, family = binomial)

s5 <- simulate(gm4,seed=101)[[1]]
s6 <- simulate(gm4,seed=101,use.u=TRUE)[[1]]

## Bernoulli
## works, but too slow
if (testLevel>2) {
    if(require("mlmRev")) {
        data(guImmun,package="mlmRev")
        g1 <- glmer(immun~kid2p+mom25p+ord+ethn+momEd+husEd+momWork+rural+pcInd81+
                    (1|comm/mom),family="binomial",data=guImmun)
        s2 <- simulate(g1)
    }
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

gm_s5 <- simulate(g1,seed=102)[[1]]

d$y <- factor(c("N","Y")[d$y+1])
g1B <- glmer(y~x+(1|f),data=d,family="binomial") ## ,tolPwrss=1e-5)
s1B <- simulate(g1B,seed=102)[[1]]
stopifnot(all.equal(gm_s5,as.numeric(s1B)-1))

## another Bernoulli
if(require("mlmRev")) {
    data(Contraception,package="mlmRev")
    gm5 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
    s3 <- simulate(gm5)
}
d$y <- rpois(nrow(d),exp(d$eta))
gm6 <- glmer(y~x+(1|f),data=d,family="poisson")
s4 <- simulate(gm6)

## simulation 'from scratch' with formulas:
## binomial
form <- formula(gm1)[-2]
gm1_s4 <- simulate(form,newdata=model.frame(gm1),
               newparams=list(theta=getME(gm1,"theta"),
               beta=fixef(gm1)),
               family=binomial,
               weights=rowSums(model.frame(gm1)[[1]]),
               seed=101)[[1]]
stopifnot(all.equal(gm1_s2,gm1_s4))


tt <- getME(gm1,"theta")
bb <- fixef(gm1)
expect_message(simulate(form,newdata=model.frame(gm1),
               newparams=list(theta=unname(tt),
               beta=fixef(gm1)),
               family=binomial,
               weights=rowSums(model.frame(gm1)[[1]]),
               seed=101),"assuming same order")
expect_error(simulate(form,newdata=model.frame(gm1),
               newparams=list(theta=setNames(tt,"abc"),
               beta=fixef(gm1)),
               family=binomial,
               weights=rowSums(model.frame(gm1)[[1]]),
               seed=101),"mismatch between")
expect_message(simulate(form,newdata=model.frame(gm1),
               newparams=list(theta=tt,
               beta=unname(bb)),
               family=binomial,
               weights=rowSums(model.frame(gm1)[[1]]),
               seed=101),"assuming same order")
expect_error(simulate(form,newdata=model.frame(gm1),
               newparams=list(theta=tt,
               beta=setNames(bb,c("abc",names(bb)[-1]))),
               family=binomial,
               weights=rowSums(model.frame(gm1)[[1]]),
               seed=101),"mismatch between")

## Gaussian
form <- formula(fm1)[-2]
s7 <- simulate(form,newdata=model.frame(fm1),
               newparams=list(theta=getME(fm1,"theta"),
               beta=fixef(fm1),
               sigma=sigma(fm1)),
               family=gaussian,
               seed=101)[[1]]
stopifnot(all.equal(s7,s1))

## TO DO: wider range of tests, including offsets ...


}
