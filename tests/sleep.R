lmod <- lFormula(Reaction ~ Days + (Days|Subject),
                 sleepstudy)
optim(c(1, 0, 1), pls, lmod = lmod,
      y = sleepstudy$Reaction)$par
mML <- lmer(Reaction ~ Days + (Days|Subject),
            sleepstudy, REML = FALSE)
getME(mML, "theta")
