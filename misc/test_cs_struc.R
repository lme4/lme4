f_install <- function(lib, branch) {
  if (dir.exists(lib)) {
    message("dir exists, skipping")
    return(invisible(NULL))
  }
  dir.create(lib)
  remotes::install_github("lme4/lme4", ref = branch, lib = lib)
}  

f_install("lib_M1", "flexSigma")
f_install("lib_M2", "flexSigma_expcs")

testfun <- function(lib) {
  library("lme4", lib = lib)
  on.exit(detach("package:lme4", unload = TRUE))
  m <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
  ret <- list(theta = getME(m, "theta"), vc = unclass(VarCorr(m)), nll = -1*c(logLik(m)))
  print(VarCorr(m))
  return(ret)
}

r1 <- testfun("lib_M1")
r2 <- testfun("lib_M2")
stopifnot(all.equal(r1$nll, r2$nll))
stopifnot(all.equal(r1$theta, r2$theta))

data("sleepstudy", package = "lme4")
gfit <- glmmTMB::glmmTMB(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)
## avoid loading glmmTMB/calling logLik method
gfit_nll <- c(gfit$obj$fn()) ## 875.9697
stopifnot(all.equal(gfit_nll, r1$nll))
gp <- gfit$fit$par
lme4:::cs_theta_to_rho(0.1629547, 2)
gp_theta <- gp[names(gp) == "theta"]
stopifnot(all.equal(unname(exp(gp_theta[1:2])),
                    unname(sqrt(diag(r1$vc$Subject))),
                    tolerance = 3e-3))

## exploration
if (FALSE) {
  devtools::load_all()
  m <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)

  trace("compute_covariance_matrix", sig = "CSCovariance", tracer = browser)
  VarCorr(m)
  object@vparameters ## [1] 0.9261518 0.2233850
  (th <- getME(m, "theta"))
  ## Subject.(Intercept) Subject.Days.(Intercept)       Subject.Days 
  ## 0.92615180               0.22338496               0.07556346
  ## not sure why corr is listed in R as ~ 0.04 rather than ~ 0.08?
  ## check to see how lambdat is actually constructed ...
  ## maybe the names are just wrong?

  obj_cs_het <- new("HeterogeneousCSCovariance", dimension = 2L)
  compute_lambdat_x(obj_cs_het, th)
  d <- 2
  sdvec <- exp(0.5*th[1:2]) ## exp(0.5*th[1:2]) ## [1] 1.588954 1.118169
  cs_theta_to_rho(th[d + 1], d)  ## 0.0377 (not 0.08)
  sdvec * sigma(m)  ## replicates the scaled SDs: [1] 40.66421 28.61597
  ## cs_theta_to_rho called with theta= 0.07556346
}

## this doesn't do what we want (treats Days as a factor in random effect?)
library(nlme)
ss <- groupedData(Reaction~1|Subject, sleepstudy)
lme(Reaction ~ Days,
    random = pdCompSymm( form = ~ Days),
    data = ss,
    method = "ML")
