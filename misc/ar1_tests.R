library(glmmTMB)
library(lme4) ## flexSigma branch!
have_asreml <- require("asreml", quietly = TRUE)
library(tinyplot)

# generate data -----------------------------------------------------------
set.seed(1)
# https://www.r-bloggers.com/2020/02/generating-correlation-matrix-for-ar1-model/
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
ndays <- 20
ngrps <- 10
y <- MASS::mvrnorm(n = 1, 
                   mu = rep(0, ndays * ngrps), 
                   kronecker(diag(ngrps), ar1_cor(ndays, 0.7))) + rnorm(ndays * ngrps)

data <- expand.grid(day = factor(1:ndays), group = LETTERS[1:ngrps]) |> 
  transform(y = y)



# fit from glmm -----------------------------------------------------------
fit_glmm <- glmmTMB(y ~  1 + ar1(0 + day|group), data, REML = TRUE)
VarCorr(fit_glmm) # pretty good given the truth! rho = 0.7, sigma = 1
#> 
#> Conditional model:
#>  Groups   Name Std.Dev. Corr       
#>  group    day1 0.90933  0.674 (ar1)
#>  Residual      1.02407


## BMB gets a different answer: bug in devel version or??
## 
## Conditional model:
##  Groups   Name Std.Dev. Corr       
##  group    day1 0.74708  0.760 (ar1)
##  Residual      1.11011

## something very 
(th <- getME(fit_glmm, "theta"))
## -0.2915797  1.1709485
c(exp(th[1]), th[2]/sqrt(1+th[2]^2))
## [1] 1.2124728 0.4780936
## ... this agrees with confint below ...

## r = t/(sqrt(1+t^2)) → r^2 = t^2/(1+t^2) → 1/r^2 = 1 + 1/t^2 →
## 1/t^2 = 1/r^2-1 → t = 1/sqrt(1/r^2-1) = r/sqrt(1-r^2)

fit_glmm2 <- update(fit_glmm, start = list(theta = c(log(1), 0.7/sqrt(1-0.7^2))))
VarCorr(fit_glmm2)

## confint(fit_glmm2)
## full object is a mess ...
confint(fit_glmm)[c("Std.Dev.day1|group", "Cor.day2.day1|group"),]
confint(fit_glmm2)[c("Std.Dev.day1|group", "Cor.day2.day1|group"),]

fit_glmm_ML <- update(fit_glmm, REML = FALSE)
VarCorr(fit_glmm_ML)
pp <- profile(fit_glmm_ML)
tinyplot(value ~ .focal, facet = ~.par,
         subset(pp, abs(.focal)<10),
         facet.args = list(free=TRUE), type = "b")
## check against CRAN version ... ??



fit_lmer <- lmer(y ~  1 + ar1(0 + day|group, hom = TRUE), data, REML = TRUE,
                 control = lmerControl(check.nobs.vs.nRE = "ignore"))
getME(fit_lmer, "theta")
logLik(fit_lmer) ## much worse?

## these don't make sense ...
##
##  group.sigma    group.rho         <NA>         <NA>         <NA>         <NA> 
##  9.864459886  1.813085901  1.295038127  2.274226690  0.485472829  0.530898250 
##         <NA>         <NA>         <NA>         <NA>         <NA>         <NA> 
##  0.506372229 -0.587073142 -0.008816628  1.131480938 -1.384474637 -1.722169315 
##         <NA>         <NA>         <NA>         <NA>         <NA>         <NA> 
## -2.133672624 -0.612020076 -0.248806760  0.310903003 -2.423565873  0.996796027 
##         <NA>         <NA> 
## -2.684895785  1.533041298 




if (have_asreml) {
  ## fit from asreml ---------------------------------------------------------
  fit_asreml <- asreml(y ~ 1, random=~ ar1(day):group, data = data, maxiter = 30)
  summary(fit_asreml)$varcomp
  ##>                   component std.error  z.ratio bound %ch
  ##> day:group         0.8304949 0.3117373 2.664085     P 0.6
  ##> day:group!day!cor 0.6730390 0.1691804 3.978232     U 0.4
  ##> units!R           1.0505278 0.2814746 3.732230     P 0.0`
}
