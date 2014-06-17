library(lme4)
m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
bruteForceHat <- function(object) {
    with(getME(object, c("Lambdat", "Lambda", "Zt", "Z", "q", "X")), {
        W <- Diagonal(x = weights(object))
        cp <- rBind(cBind(Lambdat%*%Zt%*%W%*%Z%*%Lambda +
                          diag(1,q,q),
                          Lambdat%*%Zt%*%W%*%X),
                    cBind(t(X)%*%W%*%Z%*%Lambda,
                          t(X)%*%W%*%X))
        mm <- cBind(Z%*%Lambda, X)
        mm %*% solve(as.matrix(cp)) %*% t(mm)
    })
}

H <- bruteForceHat(m)
all.equal(diag(H), unname(lme4:::hatvalues.merMod(m)))
