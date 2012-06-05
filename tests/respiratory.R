## Data originally from Davis 1991 Stat. Med., as packaged in geepack
## and transformed (center, id -> factor, idctr created, levels labeled)
library(lme4)

load("respiratory.RData")
m_glmer_4.L <- glmer(outcome~center+treat+sex+age+baseline+(1|idctr),
                     family=binomial,data=respiratory)
## works for nAGQ={2,3,5}, fails otherwise
m_glmer_4.GHQ8 <- glmer(outcome~center+treat+sex+age+baseline+(1|idctr),
                        family=binomial,data=respiratory,nAGQ=5)
