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

d <- data.frame(
  g = sample(c("A","B","C","D","E"), 250, replace=TRUE),
  y1 = runif(250, max=100),
  y2 = sample(c(0,1), 250, replace=TRUE)
)

fm1 <- lmList(y1 ~ 1 | g, data=d)
coef(fm1)
confint(fm1)

fm2 <- lmList(y2 ~ 1 | g, data=d, family=binomial)
confint(fm2)

## FIXME: methods(class="lmList") shows a bunch of methods inherited from nlme
##    that will probably fail ... is there a way to hide these/not import them?

