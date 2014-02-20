# presumably correct condVar's from lme4.0
library(lme4.0)
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
as.numeric(attr(ranef(gm,postVar=TRUE)$herd,"postVar"))
##  [1] 0.12128867 0.13363275 0.08839850 0.17337928 0.12277914 0.14436663
##  [7] 0.10658333 0.10309812 0.21289738 0.13740279 0.09555677 0.19460241
## [13] 0.14808316 0.12631006 0.15816769

fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
as.numeric(attr(ranef(fm,postVar=TRUE)$Batch,"postVar"))
## [1] 362.3583 362.3583 362.3583 362.3583 362.3583 362.3583


# but lme4 gets it wrong
library(lme4)

gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
as.numeric(attr(ranef(gm,condVar=TRUE)$herd,"postVar"))
##  [1] 0.03566755 0.04328753 0.01894730 0.07286093 0.03654278 0.05051434
##  [7] 0.02754540 0.02577625 0.10983977 0.04575869 0.02213768 0.09178386
## [13] 0.05314387 0.03868552 0.06063020

fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
as.numeric(attr(ranef(fm,condVar=TRUE)$Batch,"postVar"))
## [1] 94.55028 94.55028 94.55028 94.55028 94.55028 94.55028

# pure R version with lme4 that closely matches lme4.0 behaviour
s2.gm <- sigma(gm)^2
Lamt.gm <- getME(gm,"Lambdat")
UtU.gm <- tcrossprod(gm@pp$LamtUt)
UtUpI.gm <- UtU.gm + Diagonal(nrow(UtU.gm))
(s2.gm*tcrossprod(tcrossprod(Lamt.gm,solve(UtUpI.gm)),Lamt.gm))@x
##  [1] 0.12125894 0.13358525 0.08837934 0.17331034 0.12273769 0.14430602
##  [7] 0.10656182 0.10308298 0.21279279 0.13734534 0.09553076 0.19451821
## [13] 0.14801431 0.12628487 0.15809627
lme4:::condVar(gm)@x # with experimental function
##  [1] 0.12125894 0.13358525 0.08837934 0.17331034 0.12273769 0.14430602
##  [7] 0.10656182 0.10308298 0.21279279 0.13734534 0.09553076 0.19451821
## [13] 0.14801431 0.12628487 0.15809627

s2.fm <- sigma(fm)^2
Lamt.fm <- getME(fm,"Lambdat")
UtU.fm <- tcrossprod(fm@pp$LamtUt)
UtUpI.fm <- UtU.fm + Diagonal(nrow(UtU.fm))
(s2.fm*tcrossprod(tcrossprod(Lamt.fm,solve(UtUpI.fm)),Lamt.fm))
## 6 x 6 sparse Matrix of class "dgCMatrix"
                                                          
## [1,] 362.3113   .        .        .        .        .     
## [2,]   .      362.3113   .        .        .        .     
## [3,]   .        .      362.3113   .        .        .     
## [4,]   .        .        .      362.3113   .        .     
## [5,]   .        .        .        .      362.3113   .     
## [6,]   .        .        .        .        .      362.3113
lme4:::condVar(fm) # with experimental function
## 6 x 6 sparse Matrix of class "dgCMatrix"
                                                          
## [1,] 362.3113   .        .        .        .        .     
## [2,]   .      362.3113   .        .        .        .     
## [3,]   .        .      362.3113   .        .        .     
## [4,]   .        .        .      362.3113   .        .     
## [5,]   .        .        .        .      362.3113   .     
## [6,]   .        .        .        .        .      362.3113


# differences between lme4.0 and lme4 with more complex models.
# these differences have to be due to the conditional variance
# estimates, because all the other estimates are very similar.
# do conditional variance estimates take uncertainty in fixed
# effects into account? i didn't think so but now i'm wondering.
library(lme4)
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
attr(ranef(fm,condVar=TRUE)$Subject,"postVar")[,,1]
##           [,1]      [,2]
## [1,] 142.61353 -23.33522
## [2,] -23.33522   6.02113
lme4:::condVar(fm)[1:2, 1:2] # experimental function
## 2 x 2 sparse Matrix of class "dgCMatrix"                        
## [1,] 145.70544 -21.44448
## [2,] -21.44448   5.31228

library(lme4.0)
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
attr(ranef(fm,postVar=TRUE)$Subject,"postVar")[,,1]
##           [,1]       [,2]
## [1,] 145.71273 -21.440414
## [2,] -21.44041   5.310927

s2 <- sigma(fm)^2
Lamt <- getME(fm, "Lambdat")
L <- getME(fm, "L")
s2 * crossprod(Lamt[,1:2], solve(L, Lamt[,1:2], system = "A"))[1:2,1:2]
