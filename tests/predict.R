library(lme4)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)

## fitted values
p0 <- predict(gm1)
## fitted values, unconditional (level-0)
p1 <- predict(gm1,REform=NA)
stopifnot(length(unique(p1))==length(unique(cbpp$period)))

matplot(cbind(p0,p1),col=1:2,type="b")
newdata <- with(cbpp,expand.grid(period=unique(period),herd=unique(herd)))
## new data, all RE
p2 <- predict(gm1,newdata)
## new data, level-0
p3 <- predict(gm1,newdata,REform=NA)
## explicitly specify RE
p4 <- predict(gm1,newdata,REform=~(1|herd))
stopifnot(all.equal(p2,p4))


p5 <- predict(gm1,type="response")
stopifnot(all.equal(p5,plogis(p0)))

matplot(cbind(p2,p3),col=1:2,type="b")

## effects of new RE levels
newdata2 <- rbind(newdata,
                  data.frame(period=as.character(1:4),herd=rep("new",4)))
stopifnot(is(try(predict(gm1,newdata2),silent=TRUE),"try-error"))
p6 <- predict(gm1,newdata2,allow.new.levels=TRUE)
stopifnot(all.equal(p2,p6[1:length(p2)]))  ## original values should match
## last 4 values should match unconditional values
stopifnot(all(tail(p6,4)==predict(gm1,newdata=data.frame(period=factor(1:4)),REform=NA)))

## multi-group model
fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
## fitted values
p0 <- predict(fm1)
## fitted values, unconditional (level-0)
p1 <- predict(fm1,REform=NA)

matplot(cbind(p0,p1),col=1:2,type="b")
newdata <- with(Penicillin,expand.grid(plate=unique(plate),sample=unique(sample)))
## new data, all RE
p2 <- predict(fm1,newdata)
## new data, level-0
p3 <- predict(fm1,newdata,REform=NA)
## explicitly specify RE
p4 <- predict(fm1,newdata,REform=~(1|plate)+(~1|sample))
p4B <- predict(fm1,newdata,REform=~(1|sample)+(~1|plate))
stopifnot(all.equal(p2,p4,p4B))

p5 <- predict(fm1,newdata,REform=~(1|sample))
p6 <- predict(fm1,newdata,REform=~(1|plate))

matplot(cbind(p2,p3,p5,p6),type="b",lty=1,pch=16)

fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
p0 <- predict(fm2)
p1 <- predict(fm2,REform=NA)
## linear model, so results should be identical patterns but smaller --
##   not including intermediate days
newdata <- with(sleepstudy,expand.grid(Days=range(Days),Subject=unique(Subject)))
p2 <- predict(fm2,newdata)
p3 <- predict(fm2,newdata,REform=NA)
p3 <- predict(fm2,newdata,REform=~(0+Days|Subject))
p4 <- predict(fm2,newdata,REform=~(1|Subject))
par(mfrow=c(2,2))
library(lattice)
tmpf <- function(data,...) {
  data$Reaction <- predict(fm2,data,...)
  xyplot(Reaction~Days,group=Subject,data=data,type="l")
}
tmpf(sleepstudy)
tmpf(sleepstudy,REform=NA)
tmpf(sleepstudy,REform=~(0+Days|Subject))
tmpf(sleepstudy,REform=~(1|Subject))



