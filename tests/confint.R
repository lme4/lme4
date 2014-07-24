library("lme4")
library("testthat")

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
c0  <- confint(fm1, method="Wald")
c0B <- confint(fm1, method="Wald",parm="Days")
expect_equal(c0[2,],c0B[1,])
expect_equal(c(c0B),c(7.437592,13.496980),tolerance=1e-6)
set.seed(101)

for (bt in c("norm","basic", "perc")) {
    confint(fm1, method="boot", boot.type=bt, nsim=10)
}
for (bt in c("stud","bca","junk")) {
    expect_error(confint(fm1, method="boot", boot.type=bt, nsim=10),
                 "should be one of")
}
if((testLevel <- lme4:::testLevel()) > 1) {
    c1 <- confint(fm1,method="profile",parm=5:6)
    expect_error(confint(fm1,method="profile",parm="Days"),
                 "must be specified as an integer")
    expect_equal(c0,c1,tolerance=2e-3)  ## expect Wald and profile _reasonably_ close
    print(c1,digits=3)
    c2  <- confint(fm1,method="boot",nsim=50,parm=5:6)
    expect_error(confint(fm1,method="boot",nsim=50,parm="Days"),
                 "must be specified as an integer")
    expect_equal(c1,c2,tolerance=2e-2)
    print(c2,digits=3)
}
if (testLevel > 10) {
    print(c1B <- confint(fm1, method="profile"))
    print(c2B <- confint(fm1, method="boot"))
    expect_equal(unname(c1B), unname(c2B), tolerance=2e-2)
}

