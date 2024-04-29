library(lme4)
devtools::load_all("~/R/pkgs/lme4")
dfo <- read.csv("misc/lme4_GH691_test.csv")
table(dfo$lc)
m1 <- lmer(
            yield ~ lc +
                (1 | g1) +
                (1 | g2) +
                (lc | g3),
            data = dfo
)
ranef(m1)
predict(m1, re.form = ~(1|g1))

Dataone <- read.csv2("misc/lme4_GH691_Dataone.csv")
Dataone$dum <- dummy(Dataone$Var4,"1")
ModLMER = lmer(Var1~(1|Var3) + (0+dummy(Var4,"1")|Var5),
               Dataone,
               control=lmerControl(check.nobs.vs.nlev="ignore",
                                   check.nobs.vs.nRE="ignore"))

m2 <- update(ModLMER, . ~ (1|Var3) + (0+dum|Var5))
## attempt to predict excluding dummy term
predict(ModLMER, re.form = ~1|Var3)
## works:: dummy() is causing the problem
predict(m2, re.form = ~1|Var3)
