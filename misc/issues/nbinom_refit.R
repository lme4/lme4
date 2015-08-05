library("lme4")

set.seed(101)
dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=1:9, rep=1:15,
          KEEP.OUT.ATTRS=FALSE)
summary(mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2))))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

## mimic glmer.nb protocol
m1 <- glmer(y ~ f1*f2 + (1|g), data=dd, family=poisson)
th <- lme4:::est_theta(m1,limit=20,eps=1e-4,trace=FALSE)
m2 <- update(m1,family=negative.binomial(theta=th))
m3 <- update(m1,family=negative.binomial(theta=0.4826813))
all.equal(m1@beta,(m1B <- refit(m1))@beta)
all.equal(m2@beta,(m2B <- refit(m2))@beta)  ## fails
all.equal(m3@beta,(m3B <- refit(m3))@beta)

sessionInfo()
