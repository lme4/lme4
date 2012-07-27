library(lme4)
library(testthat)
cake2 <- rbind(cake,tail(cake,1))
cake2[nrow(cake2),"angle"] <- NA
fm0 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
fm1 <- update(fm0,data=cake2)
expect_that(update(fm1,na.action=na.fail),throws_error("missing values in object"))
fm1_omit <- update(fm1,na.action=na.omit)
## check equal:
expect_true(all.equal(fixef(fm0),fixef(fm1)))
expect_true(all.equal(VarCorr(fm0),VarCorr(fm1)))
expect_true(all.equal(ranef(fm0),ranef(fm1)))
## works, but doesn't make much sense
fm1_pass <- update(fm1,na.action=na.pass)
expect_true(all(is.na(fitted(fm1_pass))))
fm1_exclude <- update(fm1,na.action=na.exclude)

## FIXME:  fails on Windows (not on Linux!)
expect_true(is.na(tail(predict(fm1_exclude),1)))

## test predict.lm
d <- data.frame(x=1:10,y=c(rnorm(9),NA))
lm1 <- lm(y~x,data=d,na.action=na.exclude)
predict(lm1)
predict(lm1,newdata=data.frame(x=c(1:4,NA)))

## Triq examples ...
m.lmer <- lmer (angle ~ temp + (1 | recipe) + (1 | replicate), data=cake)

# Create new data frame with some NAs in fixed effect
cake2                   <- cake
cake2$temp[1:5]         <- NA

p1_pass <- predict(m.lmer, newdata=cake2, REform=NA, na.action=na.pass)
expect_true(length(p1_pass)==nrow(cake2))
expect_true(all(is.na(p1_pass[1:5])))
p1_omit <- predict(m.lmer, newdata=cake2, REform=NA, na.action=na.omit)
p1_exclude <- predict(m.lmer, newdata=cake2, REform=NA, na.action=na.exclude)
expect_true(length(p1_omit)==nrow(na.omit(cake2)))
expect_true(all.equal(p1_exclude,p1_omit))
expect_that(predict(m.lmer, newdata=cake2, REform=NA, na.action=na.fail),
            throws_error("missing values in object"))

## now try it with REform==NULL
p2_pass <- predict(m.lmer, newdata=cake2, REform=NULL, na.action=na.pass)
expect_true(length(p2_pass)==nrow(cake2))
expect_true(all(is.na(p2_pass[1:5])))
p2_omit <- predict(m.lmer, newdata=cake2, REform=NULL, na.action=na.omit)
p2_exclude <- predict(m.lmer, newdata=cake2, REform=NULL, na.action=na.exclude)
expect_true(length(p2_omit)==nrow(na.omit(cake2)))
expect_true(all.equal(p2_exclude,p2_omit))
expect_that(predict(m.lmer, newdata=cake2, REform=NULL, na.action=na.fail),
            throws_error("missing values in object"))


## experiment with NA values in random effects -- should get
## treated 
cake3 <- cake
cake3$replicate[1:5] <- NA
expect_that(predict(m.lmer, newdata=cake3, REform=NULL),
            throws_error("NAs are not allowed in prediction data"))
p4 <- predict(m.lmer, newdata=cake3, REform=NULL, allow.new.levels=TRUE)
p4B <- predict(m.lmer, newdata=cake3, REform=~1|recipe, allow.new.levels=TRUE)
expect_true(all.equal(p4[1:5],p4B[1:5]))
p4C <- predict(m.lmer, newdata=cake3, REform=NA)
