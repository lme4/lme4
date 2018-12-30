## testing whether lme4 and nlme play nicely.  Only known issue
## is lmList-masking ...
library("lme4")
library("nlme")
fm1_lmer <- lmer(Reaction ~ Days   +       (Days|Subject), sleepstudy)
fm1_lme  <- lme (Reaction ~ Days, random = ~Days|Subject,  sleepstudy)
## variance-covariance matrices: annoyingly different structures
vc_lmer <- VarCorr(fm1_lmer)
vc_lme  <- VarCorr(fm1_lme, rdig = 8)
suppressWarnings(storage.mode(vc_lme) <- "numeric")# 2 NAs
vc_lmerx <- c(diag(vc_lmer[[1]]), attr(vc_lmer[[1]],"correlation")[1,2])
vc_lmex  <- c( vc_lme[1:2,1],     vc_lme[2,3])
stopifnot(
    all.equal(vc_lmex, vc_lmerx, tolerance= 4e-4) # had 3e-5, now see 0.000296
  , ## fixed effects (much easier) :
    all.equal(fixef(fm1_lmer), fixef(fm1_lme)) # 3.6e-15
   ,
    all.equal(unname(unlist(unclass(ranef(fm1_lmer)))),
              unname(unlist(unclass(ranef(fm1_lme)))),
              tolerance = 2e-4) # had 2e-5, now see 8.41e-5
)

fm1L_lme  <- nlme::lmList(distance ~ age | Subject, Orthodont)
fm1L_lmer <- lme4::lmList(distance ~ age | Subject, Orthodont)
stopifnot(all.equal(fixef(fm1L_lmer),
                    fixef(fm1L_lme)))
sm1L_e  <- summary(fm1L_lme)
sm1L_er <- summary(fm1L_lmer)
stopifnot(
    all.equal(coef(sm1L_e),
              coef(sm1L_er), tol=1e-12)# even tol=0 works on some Lnx 64b
)

## FIXME: test opposite order
