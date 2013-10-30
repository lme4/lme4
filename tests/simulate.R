library(lme4)

(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)
L <- load(system.file("testdata/lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))
fm1 <- fit_sleepstudy_1

s1 <- simulate(fm1,seed=101)[[1]]
s2 <- simulate(fm1,seed=101,use.u=TRUE)

## binomial (2-column and prob/weights)
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_3

gm1_s1 <- simulate(gm1,seed=101)[[1]]
gm1_s2 <- simulate(gm2,seed=101)[[1]]
stopifnot(all.equal(gm1_s1[,1]/rowSums(gm1_s1),gm1_s2))
gm1_s3 <- simulate(gm1,seed=101,use.u=TRUE)

## binomial (factor): Kubovy bug report 1 Aug 2013
d <- data.frame(y=factor(rep(letters[1:2],each=100)),
                f=factor(rep(1:10,10)))
g1 <- glmer(y~(1|f),data=d,family=binomial)
invisible(simulate(g1,nsim=10))

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
gm5 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
s3 <- simulate(gm5)

d$y <- rpois(nrow(d),exp(d$eta))
gm6 <- glmer(y~x+(1|f),data=d,family="poisson")
s4 <- simulate(gm6)

## simulation 'from scratch' with formulas
form <- formula(gm1)[-2]
gm1_s4 <- simulate(form,newdata=model.frame(gm1),
               newparams=list(theta=getME(gm1,"theta"),
               beta=fixef(gm1)),
               family=binomial,
               weights=rowSums(model.frame(gm1)[[1]]),
               seed=101)[[1]]
stopifnot(all.equal(gm1_s2,gm1_s4))

## TO DO: wider range of tests, including offsets ...

