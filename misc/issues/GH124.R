require(lme4)

dds <- data.frame(f=factor(rep(1:20,each=20)))
dds$y <- suppressMessages(simulate(~1+(1|f),
                 family="gaussian",
                 newparams=list(beta=0,theta=1,sigma=1),
                 newdata=dds,seed=101)[[1]])
## create a difference, to test effects of weighting
dds$y[1:200] <- dds$y[1:200]+5
m0 <- lmer(y~1+(1|f),data=dds)
tmpf <- function(x) c(fixef(x),v=c(VarCorr(x)[[1]]),sigma=sigma(x))
bootsum <- function(b) {
    se <- sqrt(diag(var(b$t)))
    bias <- colMeans(b$t)-b$t0
    data.frame(std_bias=bias/b$t0,Z=bias/se)
}
## simple subset
tmpf(m0 <- lmer(y~1+(1|f),data=dds[1:200,]))
bootsum(b0 <- bootMer(m0,FUN=tmpf,nsim=100,seed=101))
## bias ~ -0.07 for intercept, -0.008 for sigma

tmpf(m1 <- lmer(y~1+(1|f),data=dds,weights=rep(1:0,each=200)))
bootsum(b1 <- bootMer(m1,FUN=tmpf,nsim=100,seed=101))
## both v and sigma are messed up

## the problem isn't with refit() ...
all.equal(sigma(m1),sigma(refit(m1,newresp=dds$y)))

## nor is the problem *entirely* with exactly-zero weights ...
## here the bias is only in sigma, not in v
tmpf(m2 <- lmer(y~1+(1|f),data=dds,weights=rep(c(1,1e-5),each=200)))
bootsum(b2 <- bootMer(m2,FUN=tmpf,nsim=100,seed=102))

dds$w <- rlnorm(400,meanlog=0,sdlog=1)
tmpf(m3 <- lmer(y~1+(1|f),data=dds,weights=w))
bootsum(b3 <- bootMer(m3,FUN=tmpf,nsim=100,seed=102))
## bias only in sigma here

## something funny about the way simulate() handles weights?
set.seed(101)
matplot(simulate(m0,nsim=100),type="p",pch=".",col=1,ylim=c(-10,10))
matpoints(simulate(m3,nsim=100),type="p",pch=".",col=2)
## at some point weights are getting substituted instead of sigma ???
simulate(m1)

## still can't tell what the problem is; miscalibration in sigma somewhere?
## or are the thetas getting weird because of the miscalibration in sigma?

bootMer(m1,FUN=tmpf,nsim=100,seed=101)
lmer(y~1+(1|f),data=dd,weights=rep(0:1,each=200))
## zero weights: infinite REML criterion/gradient contains NAs

lmer(y~1+(1|f),data=dd,weights=rep(c(1e-5,1),each=200))
fixef(m2 <- lmer(y~1+(1|f),data=dd,weights=rep(0:1,each=200)))
bootMer(m2,FUN=tmpf,nsim=20,seed=101)


## ORIGINAL EXAMPLE

if (FALSE) {
    ## dd <- read.csv2(file("http://topicostropicais.net/bau/data.csv"))
    download.file("http://topicostropicais.net/bau/data.csv",dest="localbaudata.csv")
    dd <- read.csv2("localbaudata.csv")
    
    m1 <- lmer(NOTA~(1|ID_ESCOLA)+(1|ID_TURMA),data=dd,weights=PESO)

    ## pre conditions seems to be OK
    ## plot(m1)
    ## qqnorm(scale(resid(m1)),ylab="Residual quantiles",col="orange")
    ## qqline(scale(resid(m1)),col="blue")  ## somewhat skewed ...

    ## reproducing error
    VarCorr(m1) # Residual = 46.5197
    sigma(m1)
    ## with "boot", sigma intervals doesn't match the estimated value 46.5197
    confint(m1,method="boot",boot.type="norm",nsim=10)
    confint(m1,method="boot",boot.type="basic")
    confint(m1,method="boot",boot.type="perc")

    ## cause seems to be 'weights'
    m2<-lmer(NOTA~(1|ID_ESCOLA)+(1|ID_TURMA),data=dd)
    VarCorr(m2)
    confint(m2,method="boot")
    confint(m2,method="profile")
}
