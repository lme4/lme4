data(Orthodont,package="nlme")
Orthodont <- as.data.frame(Orthodont)
library(lme4)
fm1 <- lmList(Reaction ~ Days | Subject, sleepstudy)
fm1 <- lmList(Reaction ~ Days | Subject, sleepstudy, pool=TRUE)
coef(fm1)
summary(fm1)
confint(fm1)
fm2 <- lmList(distance ~ age | Subject, Orthodont)
coef(fm2)

d <- data.frame(
  g = sample(c("A","B","C","D","E"), 250, replace=TRUE),
  y1 = runif(250, max=100),
  y2 = sample(c(0,1), 250, replace=TRUE)
)

fm1 <- lmList(y1 ~ 1 | g, data=d)
coef(fm1)
confint(fm1)

fm2 <- lmList(y2 ~ 1 | g, data=d, family=binomial)
confint(fm2)


fm3 <- lmList(cbind(incidence, size - incidence) ~ period|herd,
             family=binomial, data=cbpp)
coef(fm3)

## this is a slightly odd example because the residual df from
##  these fits are in fact zero ...  so pooled.SD fails, as it should

## FIXME: methods(class="lmList") shows lots of methods inherited from nlme
##    that will probably fail ...
##  hide/fix these?

## library(reshape2)
## library(ggplot2)
## ggplot(melt(as.matrix(coef(fm3))),
##       aes(value,Var2,colour=factor(Var1)))+
##    geom_point()+
##       geom_path(aes(group=factor(Var1)))
       

if (FALSE) {
    for (i in c(unclass(methods(class="lmList")))) {
        method <- gsub("\\.lmList","",i)
        cat(method,"\n")
        ## do.call(,fm3)
        ## argh; do.call("coef",fm3) and coef(fm3) behave differently
        try(eval(parse(text=paste0(method,"(fm3)"))))
    }
}    
