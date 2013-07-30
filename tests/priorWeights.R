library("lme4")

set.seed(2)
n <- 40
w <- runif(n)
x <- runif(n)
g <- factor(sample(1:10,n,replace=TRUE))
Z <- model.matrix(~g-1);
y <- Z%*%rnorm(ncol(Z)) + x + rnorm(n)/w^.5
m <- lmer(y~x+(1|g),weights=w, REML = TRUE)

fixef_lme4.0 <- c(-0.730654, 2.028954)
stopifnot(all.equal(unname(fixef(m)), fixef_lme4.0, tol = 10^-3))

sigma_lme4.0 <- 1.736143
stopifnot(all.equal(sigma(m), sigma_lme4.0, tol = 10^-3))

Sigma_lme4.0 <- 2.356705
stopifnot(all.equal(as.vector(VarCorr(m)$g), Sigma_lme4.0, tol = 10^-3))

SE_lme4.0 <- c(0.9507008, 1.3765086)
stopifnot(all.equal(as.vector(summary(m)$coefficients[,2]), SE_lme4.0, tol = 10^-3))
