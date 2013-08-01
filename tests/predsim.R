## compare range, average, etc. of simulations to
## conditional and unconditional prediction

library(lme4)
do.plot <- FALSE

fm1 <- lmer(Reaction~Days+(1|Subject),sleepstudy)
set.seed(101)
pp <- predict(fm1)
rr <- range(usim2 <- simulate(fm1,1,use.u=TRUE)[[1]])
stopifnot(all.equal(rr,c(159.3896,439.1616),tol=1e-6))
if (do.plot) {
    plot(,ylim=rr)
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
stopifnot(all.equal(unname(ss_sum[,2]),pp,tolerance=5e-3))

## population-level prediction
pp2 <- predict(fm1,REform=NA)
ss2 <- simulate(fm1,1000,use.u=FALSE)
ss_sum2 <- t(apply(ss2,1,quantile,c(0.025,0.5,0.975)))

if (do.plot) {
    plot(pp2)
    matlines(ss_sum2,col=c(1,2,1),lty=c(2,1,2))
}

stopifnot(all.equal(unname(ss_sum2[,2]),unname(pp2),tol=8e-3))
