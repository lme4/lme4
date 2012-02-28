## tests of getME: are names correct?

library(lme4Eigen)
fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
getME(fm1,"beta")
getME(fm1,"theta")

fm2 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
getME(fm2,"beta")
getME(fm2,"theta")

getME(lmer(Reaction ~ Days + (Days|Subject), sleepstudy),"theta")
getME(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
           sleepstudy),"theta")
