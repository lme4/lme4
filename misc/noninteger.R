
library(lme4)

library(glmmTMB)
m1 <- glmer(Reaction ~ Days + (1 | Subject), family = poisson, sleepstudy)
m2 <- glmmTMB(Reaction ~ Days + (1 | Subject), family = poisson, sleepstudy)

## hatvalues:
## partial(y-hat)/partial(y)
hatvalues(m1)
hatvalues(m2)
eps <- 1e-4
## ss <- transform(sleepstudy, Reaction = Reaction + c(eps, rep(0, nrow(sleepstudy)-1)))
ss <- within(sleepstudy, Reaction[1] <- Reaction[1] + eps)
yhat_old <- predict(m1)[1]
yhat_new <- predict(update(m1, data = ss))[1]
