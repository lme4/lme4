library(lme4)
fmVA <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
              (1|id) + (1|item), family = binomial, data =
              VerbAgg)
source("allFit.R")
aa.VA <- allFit(fmVA)
save("fmVA","aa.VA",file="VA.RData")

