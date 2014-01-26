## testing whether lme4 and nlme play nicely.  Only known issue
## is lmList-masking ...
library("lme4")
library("nlme")
fm1_lmer <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm1_lme <- lme(Reaction ~ Days, random=~Days|Subject, sleepstudy)
## variance-covariance matrices: annoyingly different structures
vc_lmer <- VarCorr(fm1_lmer)
vc_lmerx <- c(diag(vc_lmer[[1]]),attr(vc_lmer[[1]],"correlation")[1,2])
vc_lme <- VarCorr(fm1_lme)
suppressWarnings(storage.mode(vc_lme) <- "numeric")
vc_lmex <- c(vc_lme[1:2,1],vc_lme[2,3])
stopifnot(all.equal(vc_lmex,vc_lmerx,tolerance=3e-5))
## fixed effects (much easier)
stopifnot(all.equal(fixef(fm1_lmer),fixef(fm1_lme)))
stopifnot(all.equal(unname(unlist(unclass(ranef(fm1_lmer)))),
                    unname(unlist(unclass(ranef(fm1_lme)))),
          tolerance=2e-5))

fm1L_lme <- nlme::lmList(distance ~ age | Subject, Orthodont)
## FIXME: lmList not working yet?
## fm1L_lmer <- lme4::lmList(distance ~ age | Subject, Orthodont)

## FIXME: test opposite order
