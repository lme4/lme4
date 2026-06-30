library(lme4)
dd  <- data.frame(g = factor(rep(1:3, each = 3)))
dd$y <- simulate(~ 1 + (1|g),
              family = poisson,
              seed = 101,
              newdata = dd,
              newparams = list(beta = 0, theta = 2))[[1]]

myident <- function(x,y, tolerance = 0) {
  isTRUE(all.equal(drop(unname(x)), drop(unname(y)), check.attributes = FALSE, tolerance = tolerance))
}

## constructive::construct(dd)
## equivalent, maybe helpful in porting to Julia etc.
dd <- data.frame(
  g = factor(rep(1:3, each = 3L)),
  y = as.integer(c(0, 0, 1, 3, 5, 4, 0, 1, 0))
)

g1 <- glmer(y ~ 1 + (1|g), data = dd, family = poisson, control = glmerControl(nAGQ0initStep = FALSE))
g2 <- update(g1, control = glmerControl(compDev = FALSE, nAGQ0initStep = FALSE),
             verbose = 100)

stopifnot(identical(VarCorr(g1), VarCorr(g2)))

dfun1 <- getME(g1, "devfun")
dfun2 <- getME(g2, "devfun")
stopifnot(!identical(environment(dfun1), environment(dfun2)))

pars <- c(1.16204380240176, -0.183842216759865)
stopifnot(identical(dfun1(pars), dfun2(pars)))

stopifnot(!get("compDev", environment(dfun2)))

## compare with lme4pureR (not identical)
remotes::install_github("lme4/lme4pureR")
library(lme4pureR)

ll <- plsform(formula(g2), data = dd, family = poisson)
devf <- do.call(pirls, c(ll, list(family=poisson)))
devf(pars)
all.equal(devf(pars), dfun1(pars), tolerance = 0)
## differs by 6e-6 (could compare in more detail)

if (FALSE) {
  debug(environment(dfun2)$pwrssUpdate)
  debug(lme4:::RglmerWrkIter)
  dfun2(pars)
}

## below I step through RglmerWrkIter, checking values against my understanding at each step

pp <- environment(dfun2)$pp
resp <- environment(dfun2)$resp
eps <- .Machine$double.eps
## reset mean (start latent variables at 0)
resp$updateMu(c(0,0,0))
beta <- pars[2]
mu <- exp(beta)

## compute sqrt(working weights) by hand
stopifnot(all(resp$mu == mu))
stopifnot(all(resp$muEta() == mu))
stopifnot(all(resp$variance() == mu))
## working weights should be mu.eta/sqrt(var) = mu/sqrt(mu) = sqrt(mu)
## not exact equality here (because of order of operations)
stopifnot(all(abs(resp$sqrtWrkWt() - sqrt(mu)) < eps))
## exact equality
stopifnot(resp$sqrtWrkWt() == mu*(1/sqrt(mu)))

## step 1 of RglmerWrkIter (easy)
pp$updateXwts(resp$sqrtWrkWt())
stopifnot(all(pp$Xwts == resp$sqrtWrkWt()))

## step 2
pp$updateDecomp()
## this updates d_RZX = d_LamtUt * d_V
## and updates d_L, d_RX
## [BMB: why?? we haven't defined d_RX, d_RZX explicitly in the glmer JSS ms.,
##  but presumably these are defined as in eq. 18 of the lmer JSS paper, i.e.
##  L = upper left block of (lower-tri) Cholesky factor; RXt = lower right block;
##  and RZXt = lower left block ...]

## also sets d_ldRX2 = 2. * d_RX.matrixLLT().diagonal().array().abs().log().sum();

## wrkResp() * sqrtWrkWt();
## wrkResp = (d_eta - d_offset) + wrkResids();
## wrkResids = (d_y - d_mu) / muEta();
## at first step, d_eta - d_offset is 0
stopifnot(all((dd$y - mu) / mu  == resp$wrkResids()))

## step 3
pp$updateRes(resp$wtWrkResp())

## wtres below is the working **response** -- not the working **residuals**
## (i.e., the argument to updateRes())

## d_V  = d_Xwts.asDiagonal() * d_X;
##  d_Vtr           = d_V.adjoint() * wtres;
##  d_Utr           = d_LamtUt * wtres;
d_X <- model.matrix(~1, data = dd)
V <- diag(pp$Xwts) %*% d_X 
stopifnot(myident(t(V) %*% resp$wtWrkResp(), pp$Vtr))
Whalf <- diag(resp$sqrtWrkWt())

vals <- getME(g2, c("Lambdat", "Zt"))
## make sure Lambda has current values
stopifnot(all(vals$Lambdat@x == pars[1]))
LamtUt <- with(vals, Lambdat %*% Zt %*% Whalf)
Utr <- LamtUt %*% resp$wtWrkResp()

stopifnot(myident(Utr, pp$Utr))

ui <- c(0,0,0) ## u value before solving (in this case)

## step 4 (because uOnly is TRUE)
cval <- pp$solveU()  ## 29.28844

## BMB: it would be nice to have an explanation of these steps ... ??

## d_delb.setZero(); // in calculation of linPred delb should be zero after solveU
## d_delu          = d_Utr - d_u0;
## d_L.solveInPlace(d_delu, CHOLMOD_P);
## d_L.solveInPlace(d_delu, CHOLMOD_L);    // d_delu now contains cu
## d_CcNumer       = d_delu.squaredNorm(); // numerator of convergence criterion
## d_L.solveInPlace(d_delu, CHOLMOD_Lt);
## d_L.solveInPlace(d_delu, CHOLMOD_Pt);

LHS <- tcrossprod(LamtUt) + diag(nrow(LamtUt))
u1_alt <- solve(LHS, Utr - ui)
u1 <- solve(LHS, LamtUt %*% resp$wtWrkResp())
stopifnot(myident(Utr-ui, LamtUt %*% resp$wtWrkResp()))  ## these are in fact identical ...

## result of solving is the same (slightly larger discrepancy since C++ code is being
## fancy and updating in place, etc. etc. (permutations probably don't matter in this case,
## LamtUt is diagonal anyway ...)
stopifnot(myident(u1, pp$delu, 2.1*eps))

## numerator of convergence is criteria is squared norm of ... ???
## d_delb          = d_RX.matrixL().solve(d_Vtr - d_RZX.adjoint() * d_delu);

stopifnot(identical(pp$delu, pp$u(1)))

## b(f): d_Lambdat.adjoint() * u(f);
## linpred: d_X * beta(f) + d_Zt.adjoint() * b(f);
## this is potentially confusing because d_X * beta(f) is the _change of Xbeta from the offset_
## (not just the naive X%*%beta)

## b(f) = Lamdat*u1

bvec <- t(vals$Lambdat) %*% u1
stopifnot(all.equal(unname(drop(bvec)), pp$b(1), tolerance = .Machine$double.eps))
stopifnot(all.equal(unname(drop(with(vals, t(Zt) %*% bvec))), pp$linPred(1), tolerance = 2.1*eps))

## sets mu = offset + linPred(1)
resp$updateMu(pp$linPred(1)) ## returns updated WRSS
## Rglmerwkiter returns penalized *sum of dev resids* -- not the same as resdev ....
resp$resDev() + pp$sqrL(1)

resp$resDev()
resp$wrss()
