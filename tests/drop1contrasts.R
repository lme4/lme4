## drop1 may not work right with contrasts: make up an example something like this ...
## options(contrasts=c("contr.sum","contr.poly"))
## drop1(fecpoiss_lm3,test="Chisq",scope=.~.)

library(lme4)
oldopts <- options(contrasts=c("contr.sum","contr.poly"))
fm1 <- lmer(Reaction~Days+(Days|Subject),data=sleepstudy)
drop1(fm1,test="Chisq")
## debug(lme4:::drop1.merMod)
drop1(fm1,test="Chisq",scope=.~.)

fm0 <- lm(Reaction~Days+Subject,data=sleepstudy)
drop1(fm0,test="Chisq",scope=.~.)
options(oldopts)  ## restore original contrasts

ff <- function() {
    lmer(Reaction~Days+(Days|Subject),data=sleepstudy)
}
drop1(ff())  ## OK because sleepstudy is accessible!


