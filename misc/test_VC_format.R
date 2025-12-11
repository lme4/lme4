library(lme4)
library(glmmTMB)

fm1 <- lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
VarCorr(fm1)

sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                sleepstudy, REML = FALSE)
VarCorr(fm1.ar1)
fm1.ar1het <- lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = FALSE),
                sleepstudy, REML = FALSE)
VarCorr(fm1.ar1het)

fm2.ar1het <- glmmTMB(Reaction ~ Daysf + hetar1(0 + Daysf | Subject),
                sleepstudy, REML = FALSE)
VarCorr(fm2.ar1het)

v1 <- unclass(VarCorr(fm1.ar1het))$Subject
v2 <- unclass(VarCorr(fm2.ar1het))$cond$Subject

library(cowplot)
library(Matrix)
plot_grid(image(Matrix(v1)), image(Matrix(v2)))
