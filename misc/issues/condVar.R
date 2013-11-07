# presumably correct condVar's from lme4.0
library(lme4.0)
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
as.numeric(attr(ranef(gm,postVar=TRUE)$herd,"postVar"))

fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
as.numeric(attr(ranef(fm,postVar=TRUE)$Batch,"postVar"))

# but lme4 gets it wrong
library(lme4)
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
as.numeric(attr(ranef(gm,condVar=TRUE)$herd,"postVar"))

fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
as.numeric(attr(ranef(fm,condVar=TRUE)$Batch,"postVar"))

# pure R version with lme4 that closely matches lme4.0 behaviour
s2.gm <- sigma(gm)^2
Lamt.gm <- getME(gm,"Lambdat")
UtU.gm <- tcrossprod(gm@pp$LamtUt)
UtUpI.gm <- UtU.gm + Diagonal(nrow(UtU.gm))
(s2.gm*crossprod(tcrossprod(Lamt.gm,solve(UtUpI.gm)),Lamt.gm))@x

s2.fm <- sigma(fm)^2
Lamt.fm <- getME(fm,"Lambdat")
UtU.fm <- tcrossprod(fm@pp$LamtUt)
UtUpI.fm <- UtU.fm + Diagonal(nrow(UtU.fm))
(s2.fm*crossprod(tcrossprod(Lamt.fm,solve(UtUpI.fm)),Lamt.fm))@x
