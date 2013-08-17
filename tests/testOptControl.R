## https://github.com/lme4/lme4/issues/59
library(lme4)
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))
NMcopy <- lme4:::Nelder_Mead

cc <- capture.output(lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
                          control=
                          lmerControl("NMcopy", 
                                      optCtrl= list(iprint=20))))
## check that printing goes through step 140 twice and up to 240 once
findStep <- function(str,n) sum(grepl(paste0("^\\(NM\\) ",n,": "),cc))
stopifnot(findStep(cc,140)==2 && findStep(cc,240)==1)
