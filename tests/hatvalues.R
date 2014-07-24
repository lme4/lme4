library(lme4)
m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
bruteForceHat <- function(object) {
    with(getME(object, c("Lambdat", "Lambda", "Zt", "Z", "q", "X")), {
        ## cp:= the cross product block matrix in (17) and (18):
        W <- Diagonal(x = weights(object))
        I <- Diagonal(q)
        A.21 <- t(X) %*% W %*% Z %*% Lambda
        cp <- rBind(cBind(Lambdat %*% Zt %*% W %*% Z %*% Lambda + I, t(A.21)),
                    cBind(A.21, t(X) %*% W %*% X))
        mm <- cBind(Z %*% Lambda, X)
        ## a bit efficient: both cp and mm are typically quite sparse
        ## mm %*% solve(as.matrix(cp)) %*% t(mm)
        mm %*% solve(cp, t(mm), sparse=FALSE)
    })
}

str(H <- bruteForceHat(m))

set.seed(7)
ii <- sample(nrow(sleepstudy), 500, replace=TRUE)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy[ii, ])

stopifnot(all.equal(diag(H),
                    unname(lme4:::hatvalues.merMod(m)), tol= 1e-14),
          all.equal(diag(bruteForceHat(m2)),
                    unname(lme4:::hatvalues.merMod(m2)), tol= 1e-14)
          )
