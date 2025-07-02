## DON'T load lme4; test is to see if glmer.nb works when
## lme4 is not loaded

set.seed(101)
dd <- data.frame(x=runif(200), f= rep(1:20, each=10))
b <- rnorm(20)
dd <- transform(dd, y = rnbinom(200, mu  = exp(1 + 2*x + b[f]), size = 2))
g <- lme4::glmer.nb(y~x + (1|f), data = dd)
stopifnot(inherits(g, "glmerMod"))
