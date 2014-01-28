library(lme4)

# Dyestuff consistent with lme4.0
lme4.0condVarDyestuff <- c(362.3583, 362.3583, 362.3583, 362.3583, 362.3583, 362.3583)
fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
lme4condVarDyestuff <- as.numeric(attr(ranef(fm,condVar=TRUE)$Batch,"postVar"))
stopifnot(all.equal(lme4.0condVarDyestuff, lme4condVarDyestuff, tolerance = 1e-3))

# sleepstudy consistent with lme4.0
lme4.0condVarsleepstudy <- matrix(c(145.71273, -21.440414,
                                    -21.44041,   5.310927), 2, 2)
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
lme4condVarsleepstudy <- attr(ranef(fm,condVar=TRUE)$Subject,"postVar")[,,1]
stopifnot(all.equal(lme4.0condVarsleepstudy, lme4condVarsleepstudy, tolerance = 1e-4))

# cbpp consistent with lme4.0
lme4.0condVarcbpp <- c(0.12128867, 0.13363275, 0.08839850, 0.17337928, 0.12277914, 0.14436663,
                       0.10658333, 0.10309812, 0.21289738, 0.13740279, 0.09555677, 0.19460241,
                       0.14808316, 0.12631006, 0.15816769)
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
lme4condVarcbpp <- as.numeric(attr(ranef(gm,condVar=TRUE)$herd,"postVar"))
stopifnot(all.equal(lme4.0condVarcbpp, lme4condVarcbpp, tolerance = 1e-3))

# return warning when multiple factor terms per factor
library(testthat)
set.seed(1)
y <- rnorm(10)
x <- rnorm(10)
g <- rep(c("A","B"),5)
m1 <- lmer(y ~ 1 + (x | g))
m2 <- lmer(y ~ 1 + (1 | g) + (0 + x | g))
ranef(m1, condVar = TRUE) # no warnings expected
expect_warning(ranef(m2, condVar = TRUE))
