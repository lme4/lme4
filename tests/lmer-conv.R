### lmer() convergence testing / monitoring / ...
##  ------------------
### The output of tests here are *not* 'diff'ed  (<==> no *.Rout.save file)
library(lme4)
## Platform - and other such info -- so we find it in old saved outputs
SysI <- Sys.info()
structure(Sys.info()[c(4,5,1:3)], class="simple.list")
sessionInfo()
## and for even more details:
c(Matrix = packageDescription("Matrix")$Built,
  lme4   = packageDescription("lme4")$Built)
if(SysI[["sysname"]] == "Linux" && require("sfsmisc")) local({
    nn <- names(.Sc <- sfsmisc::Sys.cpuinfo())
    nn <- names(.Sc <- .Sc[nn != "flags"])
    print(.Sc[grep("\\.[0-9]$", nn, invert=TRUE)])
})

## convergence on boundary warnings
load(system.file("external/test3comp.rda", package = "Matrix"))
b3 <- lmer(Y3 ~ (1|Sample) + (1|Operator/Run), test3comp, verb = TRUE)

if (isTRUE(try(data(Early, package = 'mlmRev')) == 'Early')) {
    Early$tos <- Early$age - 0.5        # time on study
    b1 <- lmer(cog ~ tos + trt:tos + (tos|id), Early, verb = TRUE)
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
