## n.b. the tests run here were not really all done in the same
## run -- the runs were done separately and stitched together --
## so it's possible there are some inconsistencies in naming etc.
## across runs with different optimizers (hopefully not)
library("lme4")
library("gridExtra")

batchfn <- "penicillin_conv.RData"
nsim <- 20
lsizevec <- seq(2,6,by=0.5)

## testing
## nsim <- 3
## lsizevec <- seq(2,4,by=0.5)

source("conv_simfuns.R")

## get dimensions/names of sim results, set up 
s0 <- simfun(size=1000,seed=1001)
res0 <- expand.grid(log10size=lsizevec,rep=seq(nsim))
resB <- setNames(as.data.frame(matrix(NA,nrow=nrow(res0),ncol=length(s0))),
                 names(s0))
res <- res0 <- cbind(res0,resB)

## run nlminb
k <- 1
for (i in seq(nsim)) {
    for (j in seq_along(lsizevec)) {
        cat(i,lsizevec[j],"\n")
        res[k,-(1:2)] <- simfun(size=round(10^lsizevec[j]),seed=1000+i,
                                control=lmerControl(optimizer="nlminbwrap"))
        k <- k+1
        save("res",file=batchfn)
    }
}
res1 <- data.frame(data="Penicillin",optimizer="nlminb",res)

## reset, run bobyqa (default)
res <- res0
k <- 1
for (i in seq(nsim)) {
    for (j in seq_along(lsizevec)) {
        cat(i,lsizevec[j],"\n")
        res[k,-(1:2)] <- simfun(size=round(10^lsizevec[j]),seed=1000+i)
        k <- k+1
        save("res",file=batchfn)
    }
}
res2 <- data.frame(data="Penicillin",optimizer="bobyqa",res)
res <- rbind(res1,res2)
save("res",file=batchfn)

batchfn <- "penicillin_conv_nloptr.RData"
## reset, run bobyqa (default)
res3 <- res0
k <- 1
for (i in seq(nsim)) {
    for (j in seq_along(lsizevec)) {
        cat(i,lsizevec[j],"\n")
        res3[k,-(1:2)] <- simfun(size=round(10^lsizevec[j]),seed=1000+i,
                                 control=lmerControl(optimizer="nloptwrap"))
        k <- k+1
        save("res3",file=batchfn)
    }
}
res3 <- data.frame(data="Penicillin",optimizer="nloptr_bobyqa",res3)
res <- rbind(res,res3)
save("res",file=batchfn)
