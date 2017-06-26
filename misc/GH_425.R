library(lme4); library(optimx)
library(stringi)
library(data.table)

set.seed(1423L)
# highly imbalanced outcome variable
y <- sample.int(2L, size= 910000, replace=T, prob= c(0.98, 0.02)) - 1L
# product biases
prod <- sample(letters, size= 910000, replace=T)
#  user biases
my_grps <- stringi::stri_rand_strings(n= 35000, length= 10)
grps <- rep(my_grps, each= 26)
x1 <- sample.int(2L, size= 910000, replace=T, prob= c(0.9, 0.1)) - 1L
x2 <- sample.int(2L, size= 910000, replace=T, prob= c(0.9, 0.1)) - 1L
x3 <- sample.int(2L, size= 910000, replace=T, prob= c(0.9, 0.1)) - 1L
x4 <- sample(LETTERS[1:5], size= 91000, replace=T)

dt <- data.table(y= y,
             prod= prod, grps= grps,
             x1= x1, x2= x2, x3= x3, x4= x4)

t1 <- system.time(lmer1 <- glmer(y ~ -1 + prod + x1 + x2 + x3 + x4 + (1|grps),
               data= dt, family= binomial(link= "logit"),
               control = glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),verbose=1)

save.image("GH_425.rda")

t2 <- system.time(lmer2 <- update(lmer1,
                                  control=glmerControl(optimizer="nloptwrap")))

save.image("GH_425.rda")
