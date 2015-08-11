library("lme4")

set.seed(101)
dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=1:9, rep=1:15,
          KEEP.OUT.ATTRS=FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

## mimic glmer.nb protocol

## basic Poisson fit
m.base <- glmer(y ~ f1*f2 + (1|g), data=dd, family=poisson)
all.equal(m.base@beta,(m.base.r <- refit(m.base))@beta)  ## OK

th <- lme4:::est_theta(m.base,limit=20,eps=1e-4,trace=FALSE)
th0 <- structure(0.482681268108477, SE = 0.0244825021248148)
th1 <- structure(0.482681277470945)
th2 <- 0.482681268108477

## NB update with raw number
m.numth1 <- update(m.base,family=negative.binomial(theta=0.4826813))
all.equal(m.numth1@beta,(m.numth1.r <- refit(m.numth1))@beta)  ## OK

## strip NB value
m.symth4 <- update(m.base,family=negative.binomial(theta=c(th)))
all.equal(m.symth4@beta,(m.symth4.r <- refit(m.symth4))@beta)  ## fails

## standard NB update with computed theta from est_theta (incl SE attribute)
m.symth <- update(m.base,family=negative.binomial(theta=th))
all.equal(m.symth@beta,(m.symth.r <- refit(m.symth))@beta)  ## not OK

## NB update with equivalent value
m.symth2 <- update(m.base,family=negative.binomial(theta=th0))
all.equal(m.symth2@beta,(m.symth2.r <- refit(m.symth2))@beta)  ##

## NB update with theta value (stored as variable, no SE) only
##  (off by -1e-8)
m.symth3 <- update(m.base,family=negative.binomial(theta=th1))
all.equal(m.symth3@beta,(m.symth3.r <- refit(m.symth3))@beta)  ## works

## strip NB value (off by 5e-16)
m.symth5 <- update(m.base,family=negative.binomial(theta=th2))
all.equal(m.symth5@beta,(m.symth5.r <- refit(m.symth5))@beta)  ## fails

## Surprised by difference between specifying as variable vs.
##  number; could it be a deep-copying problem?


## try with different optimizer (red herring??)
## 
##  m4 <- update(m1,family=negative.binomial(theta=th),
##              control=glmerControl(optimizer="bobyqa"))
## all.equal(m4@beta,(m4B <- refit(m4))@beta)

refit.all <- ls(pattern="^m\\..*\\.r$")
refit.all <- mget(refit.all)
(rtab <- do.call(cbind,lapply(refit.all,fixef)))

sessionInfo()

