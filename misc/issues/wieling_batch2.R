source("wieling_batchfuns.R")

load("dialectNL2.rda")

sampfrac <- 0.1
ntot <- nrow(dialectNL)

set.seed(101)
## random subset of 10% :
subdat <- dialectNL[sample(round(sampfrac*ntot)),]

fitList0 <- fitLme4.0()
fitList1 <- lapply(argList, do.call, what=fitLme4)
save("fitList0","fitList1",
     file=batchfn)

