## test coefficient extraction in the case where RE contain
## terms that are missing from the FE ...
set.seed(101)
d <- data.frame(resp=runif(100),
                var1=factor(sample(1:5,size=100,replace=TRUE)),
                var2=runif(100),
                var3=factor(sample(1:5,size=100,replace=TRUE)))
library(lme4)
mix1 <- lmer(resp ~ 0 + var1 + var1:var2 + (1|var3), data=d)
coef(mix1)
