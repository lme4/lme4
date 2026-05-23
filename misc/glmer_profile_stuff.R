library(lme4)
cbpp2 <- read.csv(system.file("vignette_data", "cbpp2.csv",
                              package = "lme4"))
cbpp2 <- transform(
  cbpp2,
  period = factor(period),
  treatment = factor(
    treatment,
    levels = c("Partial/null", "Complete", "Unknown")
  ),
  obs = factor(seq(nrow(cbpp2))),
  avg_size = avg_size / sd(avg_size, na.rm = TRUE)
)
form <- incidence / size ~ period + treatment + avg_size + (1 | herd)        
m1 <- glmer(form, family = binomial, data = cbpp2, weights = size)

p_fun <- function(x, devtol = 1e-2, ...) {
  cat(deparse(getCall(x)$formula), "\n")
  profile(x, signames = FALSE,
          devtol = devtol,
          maxpts = 250,
          parallel = "multicore", ncpus = 10, ...)
}

aa <- allFit(m1) ## can we get close enough to sval?
## no, all values are ending up at 181.8087
with(summary(aa), sort(-2*llik))

## all we need is this?
m3 <- update(m1, control = glmerControl(nAGQ0initStep = FALSE))
p4 <- p_fun(m3)
confint(p4)

p0 <- p_fun(m1, verbose = TRUE)
p1 <- p_fun(m1, verbose = TRUE, devtol = 1e-3)

##pp <- p_fun(m1, verbose = TRUE)
## confint(pp)

## should complain about parameter out of bounds!!
p2 <- p_fun(m1, verbose = TRUE, which = 8)

sval <- list(theta = 0.558322689720326,
             beta = c(-1.00490338403396,-0.986625112550775,
                      -1.12553227684992,-1.56156926370503,
                      -0.375885923006216,-0.682742580643991,-0.0365981658810725))
m2 <- update(m1, start = sval, control = glmerControl(nAGQ0initStep = FALSE))

## old deviance  181.8087 ,
##  new deviance  181.8084 ,

-2*c(logLik(m1))
-2*c(logLik(m2))
p2 <- p_fun(m2, which = 8, verbose = TRUE)
confint(p2)
dd <- as.data.frame(p2)
plot(.zeta^2 ~ .focal, type = "b", data = dd)
devtools::load_all()
debug(profile.merMod)

## devfun2 is  a function that computes -2*nll as a function of std devs, correlations, fixed effect parameters
## in this example (GLMM with only scalar REs) devfun and devfun2 should behave the same
