library(lme4)
library(lattice)
library(testthat)
do.plots <- FALSE

L <- load(system.file("testdata/lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

if (getRversion()>"3.0.0") {
    ## saved fits are not safe with old R versions

gm1 <- fit_cbpp_1
## glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##             data = cbpp, family = binomial)

## fitted values
p0 <- predict(gm1)
p0B <- predict(gm1,newdata=cbpp)
expect_true(all.equal(p0,p0B,tolerance=2e-5)) ## FIXME: why not closer?
## fitted values, unconditional (level-0)
p1 <- predict(gm1,ReForm=NA)
expect_true(length(unique(p1))==length(unique(cbpp$period)))


if (do.plots) matplot(cbind(p0,p1),col=1:2,type="b")
newdata <- with(cbpp,expand.grid(period=unique(period),herd=unique(herd)))
## new data, all RE
p2 <- predict(gm1,newdata)
## new data, level-0
p3 <- predict(gm1,newdata,ReForm=NA)
## explicitly specify RE
p4 <- predict(gm1,newdata,ReForm=~(1|herd))
expect_true(all.equal(p2,p4))


p5 <- predict(gm1,type="response")
expect_true(all.equal(p5,plogis(p0)))

if (do.plots) matplot(cbind(p2,p3),col=1:2,type="b")

## effects of new RE levels
newdata2 <- rbind(newdata,
                  data.frame(period=as.character(1:4),herd=rep("new",4)))
expect_true(is(try(predict(gm1,newdata2),silent=TRUE),"try-error"))
p6 <- predict(gm1,newdata2,allow.new.levels=TRUE)
expect_true(all.equal(p2,p6[1:length(p2)]))  ## original values should match
## last 4 values should match unconditional values
expect_true(all(tail(p6,4)==predict(gm1,newdata=data.frame(period=factor(1:4)),ReForm=NA)))

## multi-group model
fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
## fitted values
p0 <- predict(fm1)
## fitted values, unconditional (level-0)
p1 <- predict(fm1,ReForm=NA)

if (do.plots) matplot(cbind(p0,p1),col=1:2,type="b")
newdata <- with(Penicillin,expand.grid(plate=unique(plate),sample=unique(sample)))
## new data, all RE
p2 <- predict(fm1,newdata)
## new data, level-0
p3 <- predict(fm1,newdata,ReForm=NA)
## explicitly specify RE
p4 <- predict(fm1,newdata,ReForm=~(1|plate)+(~1|sample))
p4B <- predict(fm1,newdata,ReForm=~(1|sample)+(~1|plate))
expect_true(all.equal(p2,p4,p4B))

p5 <- predict(fm1,newdata,ReForm=~(1|sample))
p6 <- predict(fm1,newdata,ReForm=~(1|plate))

if (do.plots) matplot(cbind(p2,p3,p5,p6),type="b",lty=1,pch=16)

fm2 <- fit_sleepstudy_2
p0 <- predict(fm2)
p1 <- predict(fm2,ReForm=NA)
## linear model, so results should be identical patterns but smaller --
##   not including intermediate days
newdata <- with(sleepstudy,expand.grid(Days=range(Days),Subject=unique(Subject)))
newdata$p2 <- predict(fm2,newdata)
newdata$p3 <- predict(fm2,newdata,ReForm=NA)
newdata$p4 <- predict(fm2,newdata,ReForm=~(0+Days|Subject))
newdata$p5 <- predict(fm2,newdata,ReForm=~(1|Subject))

## reference values from an apparently-working run
refval <- structure(list(Days = c(0, 9, 0, 9, 0, 9), Subject = structure(c(1L, 
1L, 2L, 2L, 3L, 3L), .Label = c("308", "309", "310", "330", "331", 
"332", "333", "334", "335", "337", "349", "350", "351", "352", 
"369", "370", "371", "372"), class = "factor"), p2 = c(253.663652396798, 
430.66001930835, 211.006415533628, 227.634788908917, 212.444742696829, 
257.61053840953), p3 = c(251.405104848485, 345.610678484848, 
251.405104848485, 345.610678484848, 251.405104848485, 345.610678484848
), p4 = c(251.405104848485, 428.401471760037, 251.405104848485, 
268.033478223774, 251.405104848485, 296.570900561186), p5 = c(253.663652396798, 
347.869226033161, 211.006415533628, 305.211989169991, 212.444742696829, 
306.650316333193)), .Names = c("Days", "Subject", "p2", "p3", 
"p4", "p5"), out.attrs = structure(list(dim = structure(c(2L, 
18L), .Names = c("Days", "Subject")), dimnames = structure(list(
    Days = c("Days=0", "Days=9"), Subject = c("Subject=308", 
    "Subject=309", "Subject=310", "Subject=330", "Subject=331", 
    "Subject=332", "Subject=333", "Subject=334", "Subject=335", 
    "Subject=337", "Subject=349", "Subject=350", "Subject=351", 
    "Subject=352", "Subject=369", "Subject=370", "Subject=371", 
    "Subject=372")), .Names = c("Days", "Subject"))), .Names = c("dim", 
"dimnames")), row.names = c(NA, 6L), class = "data.frame")

expect_true(all.equal(head(newdata),refval))

library(lattice)
tmpf <- function(data,...) {
  data$Reaction <- predict(fm2,data,...)
  if (do.plots) xyplot(Reaction~Days,group=Subject,data=data,type="l")
}
tmpf(sleepstudy)
tmpf(sleepstudy,ReForm=NA)
tmpf(sleepstudy,ReForm=~(0+Days|Subject))
tmpf(sleepstudy,ReForm=~(1|Subject))

## from 'Colonel Triq': examples using *fewer* random effect levels
##  than in original data set
m <- lmer(angle ~ temp + recipe + (1 | replicate), data=cake)
summary(m)

# replicate 1 only appears in rows 1:18.
rownames(cake[cake$replicate==1,])

predict(m, newdata=cake[-1:-17,], ReForm=~ (1 | replicate))
predict(m, newdata=cake[-1:-18,], ReForm=NA)
predict(m, newdata=cake[-1:-18,], ReForm=~ (1 | replicate)) 
predict(m, newdata=cake[-1:-18,], ReForm=~ (1 | replicate), allow.new.levels=TRUE)


}
