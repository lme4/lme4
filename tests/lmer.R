## Check that quasi families throw an error
try(gm1 <- lmer(cbind(incidence, size - incidence) ~ period + (1|herd),
                data = cbpp, family = quasibinomial))
try(gm1 <- lmer(incidence ~ period + (1|herd),
                data = cbpp, family = quasipoisson))
try(gm1 <- lmer(incidence ~ period + (1|herd),
                data = cbpp, family = quasi))

# check bug from Kevin Buhner
X <- data.frame(y=runif(10), x=rnorm(10), z=sample(c("A","B"), 10, TRUE))
lmer(log(y) ~ x | z, data=X)

