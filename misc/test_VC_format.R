library(lme4)
library(glmmTMB)

fm1 <- lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
VarCorr(fm1)

sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                sleepstudy, REML = FALSE)

VarCorr(fm1.ar1)

## FIXME: this should print all std devs/variances (should internally call the print class for hetar1 (VarCorr should internally return
## an object of class "vcmat_hetar1" rather than "vcmat_ar1"
## what is get_sd doing and how are the glmmTMB structures different/what needs to be
## done differently?

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


fm1.cs <- lmer(Reaction ~ Daysf + cs(0 + Daysf | Subject, hom = TRUE), 
               sleepstudy, REML = FALSE)
fm2.cs <- glmmTMB(Reaction ~ Daysf + homcs(0 + Daysf | Subject),
               sleepstudy, REML = FALSE)

VarCorr(fm1.cs)
VarCorr(fm2.cs)
## FIXME: lmer is printing all of the std devs even though we specified hom = TRUE ....
## FIXME: glmmTMB might not have a method for homcs, so defaulting to printing everything ... ugh ...
