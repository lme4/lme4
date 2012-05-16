## data set and formula extracted from ?prLogisticDelta example
##   (Thailand, clustered-data) in prLogistic package
load("prLogistic.RData")
library(lme4)
glmer(rgi~  sex + pped + (1|schoolid), 
                data = dataset, family=binomial)

