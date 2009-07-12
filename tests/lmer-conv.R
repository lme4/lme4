### lmer() convergence testing / monitoring / ...
##  ------------------
### The output of tests here are *not* 'diff'ed  (<==> no *.Rout.save file)
library(lme4)

## convergence on boundary warnings
load(system.file("external/test3comp.rda", package = "Matrix"))
b3 <- lmer(Y3 ~ (1|Sample) + (1|Operator/Run), test3comp, verb = TRUE)

if (isTRUE(try(data(Early, package = 'mlmRev')) == 'Early')) {
    Early$tos <- Early$age - 0.5        # time on study
    b1 <- lmer(cog ~ tos + trt:tos + (tos|id), Early, verb = TRUE)
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
