library(testthat)
## require(agridat)
## dat <- gotway.hessianfly

## don't actually use gotway_hessianfly_fit or gotway_hessianfly_prof,
## so we should be OK even with R< 3.0.1
load(system.file("testdata","gotway_hessianfly.rda",package="lme4"))
# Block random.  See Glimmix manual, output 1.18.
# Note: (Different parameterization)

## require("lme4.0")
## fit2 <- glmer(cbind(y, n-y) ~ gen + (1|block), data=dat, family=binomial)
## params <- list(fixef=fixef(fit2),theta=getME(fit2,"theta"))
## detach("package:lme4.0")
lme4.0fit <- structure(list(fixef = structure(c(1.50345713031203, -0.193853259383803, 
-0.540808391060274, -1.43419379979154, -0.203701042949808, -0.978322555343941, 
-0.604078624475678, -1.67742449813309, -1.39842466673692, -0.681709344788684, 
-1.46295367186169, -1.45908310198959, -3.55285756517073, -2.50731975980307, 
-2.08716296677356, -2.96974270029992), .Names = c("(Intercept)", 
"genG02", "genG03", "genG04", "genG05", "genG06", "genG07", "genG08", 
"genG09", "genG10", "genG11", "genG12", "genG13", "genG14", "genG15", 
"genG16")), theta = structure(0.0319087494293615, .Names = "block.(Intercept)")), .Names = c("fixef", 
"theta"))

## start doesn't work because we don't get there
library(lme4)
m1 <- glmer(cbind(y, n-y) ~ gen + (1|block), data=gotway.hessianfly,
            family=binomial)
m1B <- update(m1,control=glmerControl(optimizer="bobyqa"))
max(abs(m1@optinfo$derivs$gradient))  ## 0.0012
max(abs(m1B@optinfo$derivs$gradient)) ## 2.03e-5
abs(m1@optinfo$derivs$gradient)/abs(m1B@optinfo$derivs$gradient)
## bobyqa gets gradients *at least* 1.64* lower

lme4fit <- list(fixef=fixef(m1),theta=getME(m1,"theta"))

## hack around slight naming differences
lme4fit$theta <- unname(lme4fit$theta)
lme4.0fit$theta <- unname(lme4.0fit$theta)
expect_equal(lme4fit,lme4.0fit,tolerance=3e-4)

## Fun stuff: visualize and alternative model

## library(ggplot2)
## dat$prop <- dat$y/dat$n
## theme_set(theme_bw())
## ggplot(dat,aes(x=gen,y=prop,colour=block))+geom_point(aes(size=n))+
##     geom_line(aes(group=block,colour=block))+
##     geom_smooth(family=binomial,aes(weight=n,colour=block,group=block),method="glm",
##                 alpha=0.1)

## dat$obs <- factor(seq(nrow(dat)))
## m2 <- glmer(cbind(y, n-y) ~ block+ (1|gen) + (1|obs), data=dat, family=binomial)

