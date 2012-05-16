data(Orthodont,package="nlme")
Orthodont <- as.data.frame(Orthodont)
library(lme4)
fm1 <- lmList(Reaction ~ Days | Subject, sleepstudy)
fm1 <- lmList(Reaction ~ Days | Subject, sleepstudy, pool=TRUE)
coef(fm1)
summary(fm1)
confint(fm1)
fm2 <- lmList(distance ~ age | Subject, Orthodont)
coef(fm2)
## FIXME: add glm example?
## FIXME: methods(class="lmList") shows a bunch of methods inherited from nlme
##    that will probably fail ... is there a way to hide these/not import them?
