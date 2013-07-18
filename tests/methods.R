library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
stopifnot(colnames(model.frame(fm1))==c("Reaction","Days","Subject"))
stopifnot(colnames(model.frame(fm1,fixed.only=TRUE))==c("Reaction","Days"))
