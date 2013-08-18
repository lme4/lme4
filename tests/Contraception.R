options(show.signif.stars = FALSE)
library(lme4)
library(nloptwrap)
#library(mlmRev)
data(Contraception, package="mlmRev")
Contraception <- within(Contraception, ch <- factor(as.numeric(as.integer(livch)>1L)))
form <- use ~ age + I(age^2) + ch + urban + (1|district)
gm0 <- glm(nobars(form), binomial, Contraception)
print(summary(gm0), corr=FALSE)
## Call: glm(formula = nobars(form), family = binomial, data = Contraception)

## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -1.454  -1.037  -0.668   1.238   1.979  

## Coefficients:
##              Estimate Std. Error z value Pr(>|z|)
## (Intercept) -0.942935   0.149730   -6.30  3.0e-10
## age          0.004967   0.007589    0.65     0.51
## I(age^2)    -0.004328   0.000692   -6.26  4.0e-10
## ch1          0.806217   0.142262    5.67  1.5e-08
## urbanY       0.766362   0.105946    7.23  4.7e-13

## (Dispersion parameter for binomial family taken to be 1)

##     Null deviance: 2590.9  on 1933  degrees of freedom
## Residual deviance: 2417.9  on 1929  degrees of freedom
## AIC: 2428

## Number of Fisher Scoring iterations: 4

dput(deviance(gm0)) # 2417.86412239273
logLik(gm0) # -1208.93206119637
beta0 <- coef(gm0)

# gm1 <- glmer(form, Contraception, binomial, nAGQ=1L)
## Generalized linear mixed model fit by maximum likelihood ['summary.merMod']
##  Family: binomial ( logit )
## Formula: use ~ age + I(age^2) + ch + urban + (1 | district) 
##    Data: Contraception 

##      AIC      BIC   logLik deviance 
##  2385.19  2418.59 -1186.59  2373.19 

## Random effects:
##  Groups   Name        Variance Std.Dev.
##  district (Intercept) 0.225    0.474   
## Number of obs: 1934, groups: district, 60

## Fixed effects:
##              Estimate Std. Error z value Pr(>|z|)
## (Intercept) -1.006374   0.167892   -5.99  2.0e-09
## age          0.006256   0.007840    0.80     0.42
## I(age^2)    -0.004635   0.000716   -6.47  9.7e-11
## ch1          0.860382   0.147353    5.84  5.3e-09
## urbanY       0.692922   0.119668    5.79  7.0e-09

glmod <- glFormula(form, Contraception, binomial)
devf <- pirls(glmod, gm0$y, gm0$linear.predictor)
devf(c(1,beta0))
beta1 <- c(-1.00637368953914, 0.00625611018864605, -0.00463527887823825, 0.860382098849009, 0.692921941313297)
theta1 <- 0.473985228740739
devf(c(theta1,beta1))
# deviance(gm1)  # 2373.18581907321
# head(getME(gm1,"u")) # c(-1.56655252109295, -0.0451095377205166,
#    0.453409972754753, 0.35265189198274, 0.252723822217955, -0.504481434155787)
# head(getME(gm1, "mu")) # c(0.16025103705494, 0.225473461055359,
#    0.451107632834752, 0.383912222123962, 0.119942471057264, 0.148334820475519)

form1 <- use ~ age + I(age^2) + ch + urban + (1|district) + (1|district:urban)

## > print(summary(gm2 <- glmer(form1,Contraception,binomial)), corr=FALSE)
## Generalized linear mixed model fit by maximum likelihood ['summary.merMod']
##  Family: binomial ( logit )
## Formula: use ~ age + I(age^2) + ch + urban + (1 | district) + (1 | district:urban) 
##    Data: Contraception 

##      AIC      BIC   logLik deviance 
##  2375.88  2414.85 -1180.94  2361.88 

## Random effects:
##  Groups         Name        Variance Std.Dev.
##  district:urban (Intercept) 0.31863  0.5645  
##  district       (Intercept) 0.00715  0.0845  
## Number of obs: 1934, groups: district:urban, 102; district, 60

## Fixed effects:
##              Estimate Std. Error z value Pr(>|z|)
## (Intercept) -1.032745   0.175593   -5.88  4.1e-09
## age          0.005876   0.007935    0.74     0.46
## I(age^2)    -0.004537   0.000724   -6.26  3.8e-10
## ch1          0.872691   0.149028    5.86  4.7e-09
## urbanY       0.764670   0.169725    4.51  6.6e-06

beta2 <- c(-1.03274535458238, 0.00587602227525672, -0.0045372320199417,
           0.872690858261632, 0.764670407439094)
theta2 <- c(0.564472021994911, 0.0845334785825269)

## dput(head(unname(getME(gm2,"u"))))
## c(-1.60622219325897, -0.981855964992823, -0.0339024704692014, 
##    0.505250674118111, -0.428344594466731, 1.10122436874433)

## dput(head(unname(getME(gm2,"mu"))))
## c(0.195157958203462, 0.26347290296241, 0.504168959587931, 0.436350391068322, 
## 0.145679693412778, 0.178092903689705)

form2 <- use ~ age + I(age^2) + ch + urban + (1 | district:urban)

## print(summary(gm3 <- glmer(form2,Contraception,binomial)), corr=FALSE)
## Generalized linear mixed model fit by maximum likelihood ['summary.merMod']
##  Family: binomial ( logit )
## Formula: use ~ age + I(age^2) + ch + urban + (1 | district:urban) 
##    Data: Contraception 

##      AIC      BIC   logLik deviance 
##  2373.88  2407.29 -1180.94  2361.88 

## Random effects:
##  Groups         Name        Variance Std.Dev.
##  district:urban (Intercept) 0.327    0.572   
## Number of obs: 1934, groups: district:urban, 102

## Fixed effects:
##              Estimate Std. Error z value Pr(>|z|)
## (Intercept) -1.032829   0.175651   -5.88  4.1e-09
## age          0.005862   0.007935    0.74     0.46
## I(age^2)    -0.004535   0.000724   -6.26  3.9e-10
## ch1          0.872716   0.149034    5.86  4.7e-09
## urbanY       0.766667   0.170601    4.49  7.0e-06

beta3 <- c(-1.03282920682185, 0.00586163322902896, -0.00453480196155543, 
           0.872715839072072, 0.766667032721774)
theta3 <- 0.571568547571509

## dput(head(unname(getME(gm3,"mu"))))
## c(0.195804532044805, 0.264187372529822, 0.505051933870551, 0.43723605274129, 
## 0.146198921135808, 0.178681337173397)

## dput(head(unname(getME(gm3,"u"))))
## c(-1.63882047908625, -1.02416872093547, -0.0343530888451908, 
## 0.510937117406129, -0.419566092717808, 1.10846888333543)
