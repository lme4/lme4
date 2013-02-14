library("lme4")
testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
c0 <- confint(fm1,method="Wald")
c0B <- confint(fm1,method="Wald",parm="Days")
print(c0,digits=3)
if (testLevel>1) {
    c1 <- confint(fm1,method="profile",parm=5:6)
    try(c1B <- confint(fm1,method="profile",parm="Days"))
    print(c1,digits=3)
    c2 <- confint(fm1,method="boot",nsim=50,parm=4:5)
    c2B <- confint(fm1,method="boot",nsim=50,parm="Days")
    print(c2,digits=3)
}
