## Example from Rune's pull request:
##>>> https://github.com/lme4/lme4/pull/144
## More "handling" in
##>>> https://github.com/lme4/lme4/issues/168

## [Rune: ]
## I have updated relevant tests and added some tests related to `anova` on
## single model fits which are affected.

## The main benefit of this change is the ability to fit models such as

## borrowed from ../../inst/tests/test-rank.R
set.seed(11)##   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d1 <- expand.grid(a=factor(1:4), b=factor(1:4), rep=1:10)
n <- nrow(d1)
d1 <- transform(d1, r = sample(1:5, size=n, replace=TRUE),
                    z = rnorm(n))
d2 <- subset(d1, !(a=="4" & b=="4")) # so that is "missing"
fm <- lmer( z ~ a*b + (1|r), data=d2)
## design is column rank deficient so dropping 1 coef
anova(fm)
## Analysis of Variance Table
##     Df Sum Sq Mean Sq F value
## a    3 3.6087 1.20290  1.3046
## b    3 4.4842 1.49473  1.6211
## a:b  8 4.6847 0.58559  0.6351

## As you see, [g]lmer signals that it drops 1 coefficient (and maybe you
## want an option to turn of that message?). The best thing would naturally
## be to simulate the `lm` behaviour and have the offending coefficient show
## up as `NA` in `print` and `summary`, but that involves much deeper changes
## to the code base:

fm
## Linear mixed model fit by REML ['lmerMod']
## Formula: z ~ a * b + (1 | r) 
##    Data: d2 
## REML criterion at convergence: 415.2903 
## Random effects:
##  Groups   Name        Std.Dev.
##  r        (Intercept) 0.06589 
##  Residual             0.98962 
## Number of obs: 150, groups: r, 5
## Fixed Effects:
## (Intercept)           a2           a3           a4           b2           b3  
##     0.07387     -0.51343      0.28435      0.32098     -0.14915     -0.10380  
##          b4        a2:b2        a3:b2        a4:b2        a2:b3        a3:b3  
##     0.25370      0.25457     -0.32327     -0.19222      0.48129     -0.49072  
##       a4:b3        a2:b4        a3:b4  
##    -0.23367      0.83704     -0.78918  

lm(z ~ a*b, data=d2)
##    ............
## Coefficients:
##    ............
##       a4:b3        a2:b4        a3:b4        a4:b4  
##    -0.24603      0.84076     -0.79106           NA  
##						  ^^^^

## Here is another 'real world' example which benefits from the change:

data(soup, package="ordinal")
soup <- within(soup,
               sureness <- as.numeric(as.character(SURENESS)))
fm2 <- lmer(sureness ~ PRODID * DAY + (1|RESP), data=soup)
## design is column rank deficient so dropping 1 coef
(am2 <- anova(fm2))
## Analysis of Variance Table
##            Df Sum Sq Mean Sq F value
## PRODID      5 638.06 127.611 41.3472
## DAY         1  19.29  19.289  6.2497
## PRODID:DAY  4  27.25   6.811  2.2070

## so now I can get the ANOVA F-test of the interaction. 
stopifnot(am2[3,"Df"] == 4,
          all.equal(am2[3,"F value"], 2.207, tol=0.01))

fm1 <- update(fm2, ~ .-PRODID:DAY)
(a12 <- anova(fm1,fm2))
## barely *not* significant :
stopifnot(a12[2,"Chi Df"] == 4,
          all.equal(a12[2,]$Pr, 0.06475849, tol = 0.01))
