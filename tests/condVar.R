lme4.0condVarDyestuff <- c(362.3583, 362.3583, 362.3583, 362.3583, 362.3583, 362.3583)
fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
lme4condVarDyestuff <- as.numeric(attr(ranef(fm,condVar=TRUE)$Batch,"postVar"))
stopifnot(all.equal(lme4.0condVarDyestuff, lme4condVarDyestuff, tol = 1e-3))

lme4.0condVarcbpp <- c(0.12128867, 0.13363275, 0.08839850, 0.17337928, 0.12277914, 0.14436663,
                       0.10658333, 0.10309812, 0.21289738, 0.13740279, 0.09555677, 0.19460241,
                       0.14808316, 0.12631006, 0.15816769)
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
lme4condVarcbpp <- as.numeric(attr(ranef(gm,condVar=TRUE)$herd,"postVar"))
stopifnot(all.equal(lme4.0condVarcbpp, lme4condVarcbpp, tol = 1e-3))
