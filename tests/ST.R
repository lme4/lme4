require(lme4)
# sorry for fitting yet another sleepstudy model in the tests
m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
ST <- getME(m, "ST")$Subject

# copied from vince dorie's simmer.R in arm:
dimension <- nrow(ST)
T <- ST
diag(T) <- rep(1, dimension)
S <- diag(diag(ST), dimension)

vc0 <- getME(m, 'sigma')^2*tcrossprod(T %*% S)
vc1 <- VarCorr(m)$Subject[,]
dimnames(vc0) <- dimnames(vc1)

all.equal(vc0, vc1, tolerance = 1e-6)
