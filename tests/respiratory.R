## Data originally from Davis 1991 Stat. Med., as packaged in geepack
## and transformed (center, id -> factor, idctr created, levels labeled)
library(lme4)

load(system.file("testdata","respiratory.RData",package="lme4"))
m_glmer_4.L <- glmer(outcome~center+treat+sex+age+baseline+(1|idctr),
                     family=binomial,data=respiratory)
## FIXME: works for nAGQ={2,3,5}, fails otherwise
m_glmer_4.GHQ8 <- glmer(outcome~center+treat+sex+age+baseline+(1|idctr),
                        family=binomial,data=respiratory,nAGQ=5)

## m_glmer_4.GHQ8 <- glmer(outcome~center+treat+sex+age+baseline+(1|idctr),
##                        family=binomial,data=respiratory,nAGQ=8,verbose=10)
