library(lme4)
library(nloptwrap)
lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
devf <- pls(lmod,sleepstudy$Reaction)
bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))[c("par","value")]
mML <- lmer(Reaction ~ Days + (Days|Subject),
            sleepstudy, REML = FALSE)
getME(mML, "theta")
deviance(mML)
