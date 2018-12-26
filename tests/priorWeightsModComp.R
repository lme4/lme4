library(lme4)
n <- nrow(sleepstudy)
options(warn = 1, # show as they happen ("false" convergence warnings)
        useFancyQuotes = FALSE)
##' remove all attributes but names
dropA <- function(x) `attributes<-`(x, list(names = names(x)))
##' transform result of "numeric" all.equal.list() to a named vector
all.eqL <- function(x1, x2, ...) {
    r <- sub("^Component ", '', all.equal(x1, x2, tolerance = 0, ...))
    r <- strsplit(sub(": Mean relative difference:", "&&", r),
                  split="&&", fixed=TRUE)
    setNames(as.numeric(vapply(r, `[`, "1.234", 2L)),
             ## drop surrounding "..."
             nm = sub('"$', '', substring(vapply(r, `[`, "nam",   1L), first=2)))
}
seedF <- function(s) {
    if(s %in% c(6, 39, 52, 57, 63, 74, 76, 86))
        switch(as.character(s)
               , "52"=, "63"=, "74" = 2
               , "6"=, "39" = 3
               , "86" =  8 # needs  4 on Lnx-64b
               , "76" = 70 # needs 42 on Lnx-64b
               , "57" = 90 # needs 52 on Lnx-64b
               )
    else if(s %in% c(1, 12, 15, 34, 36, 41, 42, 43, 49, 55, 59, 67, 80, 85)) ## seeds 41,59, .. 15
        1.0
    else ## seeds 22, 20, and better
        0.25
}
## be fast, running only 10 seeds by default:
sMax <- if(lme4:::testLevel() > 1) 99L else 9L
mySeeds <- 0L:sMax

lapply(setNames(,mySeeds), function(seed) {
    cat("\n------ random seed =", seed, "---------\n")
    set.seed(seed)
    v <- rpois(n,1) + 1
    w <- 1/v
    cat("weights w:\n")
    fm1    <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE, weights = w); cat("..2:\n")
    fm2    <- lmer(Reaction ~ Days + (1    | Subject), sleepstudy, REML=FALSE, weights = w)
    cat("weights w*10:\n")
    fm1.10 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE, weights = w*10);cat("..2:\n")
    fm2.10 <- lmer(Reaction ~ Days + (1    | Subject), sleepstudy, REML=FALSE, weights = w*10)
    ##
    ano12... <- dropA(anova(fm1,    fm2   ))
    ano12.10 <- dropA(anova(fm1.10, fm2.10))
    print(aEQ <- all.eqL(ano12..., ano12.10)) # showing differences
    if(!exists("notChisq"))
	notChisq <<-
	    local({ n <- names(ano12...)
		grep("Chisq", n, value=TRUE, fixed=TRUE, invert=TRUE) })
    stopifnot(
        all.equal(ano12...$Chisq,
                  ano12.10$Chisq, tol = 1e-6 * seedF(seed))
       ,
        all.equal(ano12...[notChisq],
                  ano12.10[notChisq], tol= 1.5e-8 * seedF(seed))
    )
    aEQ
}) -> rallEQ

cat("=====================================\n")

rallEQ <- t(simplify2array(rallEQ))
notChisq <- intersect(notChisq, colnames(rallEQ))
## sort according to "severity":
srallEQ <- rallEQ[with(as.data.frame(rallEQ), order(AIC, Chisq)), ]
round(log10(srallEQ), 2)
saveRDS(srallEQ, "priorWeightsMod_relerr.rds")

if(!dev.interactive(orNone=TRUE)) pdf("priorWeightsMod_relerr.pdf")

matplot(mySeeds, log10(srallEQ), type="l", xlab=NA) ; grid()
legend("topleft", ncol=3, bty="n",
       paste(1:6, colnames(srallEQ), sep = ": "), col=1:6, lty=1:6)
tolD <- sqrt(.Machine$double.eps) # sqrt(eps_C)
abline(h = log10(tolD), col = "forest green", lty=3)
axis(4, at=log10(tolD), label=quote(sqrt(epsilon[c])), las=1)
LRG <- which(srallEQ[,"AIC"] > tolD)
text(LRG, log10(srallEQ[LRG, "AIC"]), names(LRG), cex = .8)

## how close are we ..
str(tF <- sapply(mySeeds, seedF))
round(sort(      rallEQ[, "Chisq"] / (tF * 1e-6  ),          decreasing=TRUE), 1)
round(sort(apply(rallEQ[,notChisq] / (tF * 1.5e-8), 1, max), decreasing=TRUE), 1)
