## compare range, average, etc. of simulations to
## conditional and unconditional prediction
library(lme4)
do.plot <- FALSE

fm1 <- lmer(Reaction~Days+(1|Subject),sleepstudy)
set.seed(101)
pp <- predict(fm1)
rr <- range(usim2 <- simulate(fm1,1,use.u=TRUE)[[1]])
stopifnot(all.equal(rr,c(159.3896,439.1616),tolerance=1e-6))
if (do.plot) {
    plot(pp,ylim=rr)
    lines(sleepstudy$Reaction)
    points(simulate(fm1,1)[[1]],col=4)
    points(usim2,col=2)
}

set.seed(101)

## conditional prediction
ss <- simulate(fm1,1000,use.u=TRUE)
ss_sum <- t(apply(ss,1,quantile,c(0.025,0.5,0.975)))
plot(pp)
matlines(ss_sum,col=c(1,2,1),lty=c(2,1,2))
stopifnot(all.equal(ss_sum[,2],pp,tolerance=5e-3))

## population-level prediction
pp2 <- predict(fm1,ReForm=NA)
ss2 <- simulate(fm1,1000,use.u=FALSE)
ss_sum2 <- t(apply(ss2,1,quantile,c(0.025,0.5,0.975)))

if (do.plot) {
    plot(pp2,ylim=c(200,400))
    matlines(ss_sum2,col=c(1,2,1),lty=c(2,1,2))
}

stopifnot(all.equal(ss_sum2[,2],pp2,tolerance=8e-3))

## predict(...,newdata=...) on models with derived variables in the random effects
## e.g. (f:g, f/g)
set.seed(101)
d <- expand.grid(f=factor(letters[1:10]),g=factor(letters[1:10]),
                 rep=1:10)
d$y <- rnorm(nrow(d))
m1 <- lmer(y~(1|f:g),d)
p1A <- predict(m1)
p1B <- predict(m1,newdata=d)
stopifnot(all.equal(p1A,p1B))
m2 <- lmer(y~(1|f/g),d)
p2A <- predict(m2)
p2B <- predict(m2,newdata=d)
stopifnot(all.equal(p2A,p2B))

## with numeric grouping variables
dn <- transform(d,f=as.numeric(f),g=as.numeric(g))
m1N <- update(m1,data=dn)
p1NA <- predict(m1N)
p1NB <- predict(m1N,newdata=dn)
stopifnot(all.equal(p1NA,p1NB))

## simulate with modified parameters
set.seed(1)
s1 <- simulate(fm1)
set.seed(1)
s2 <- simulate(fm1,newdata=model.frame(fm1),
               newparams=getME(fm1,c("theta","beta","sigma")))
all.equal(s1,s2)

fm0 <- update(fm1,.~.-Days)
##
## sim() -> simulate() -> refit() -> deviance
##

## predictions and simulations with offsets

set.seed(101)
d <- data.frame(y=rpois(100,5),x=rlnorm(100,1,1),
                f=factor(sample(10,size=100,replace=TRUE)))
gm1 <- glmer(y~offset(log(x))+(1|f),data=d,
             family=poisson)
s1 <- simulate(gm1)
