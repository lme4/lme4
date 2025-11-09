
test_that("unit test for unstructured covariances", {
  
  for (nc in c(0:4, 16L, 64L, 256L)){
    
    x.us <- new("Covariance.us", nc = nc, simulate = TRUE)
    
    ## Basic sanity tests
    expect_s4_class(x.us, "Covariance.us")
    expect_true(validObject(x.us))
    
    ## par should be the same as theta
    expect_equal(getPar(x.us), getTheta(x.us))
    
    ## checking for specific par/theta lengths
    expect_equal(getParLength(x.us), (nc * (nc+1)/2))
    expect_equal(getParLength(x.us), getThetaLength(x.us))
    expect_equal(length(getThetaIndex(x.us)), 
                 getThetaIndexLength(x.us))
    expect_equal(getLambdat.dp(x.us), seq_len(nc))
    
    ## Testing getUpper() and getLower() in a different manner
    expect_equal(getUpper(x.us),
                 rep(Inf, getThetaIndexLength(x.us)))
    
    if(nc > 1L){
      test_low <- rep(-Inf, nc * (nc + 1) / 2)
      diag_positions <- c(1, cumsum(nc:2) + 1)
      test_low[diag_positions] <- 0
    } else if (nc == 1L){
      test_low <- 0
    } else if (nc == 0L){
      test_low <- numeric(0)
    }
    expect_equal(getLower(x.us), test_low)
    
    ## Testing getLambdat.i in a different manner
    expect_equal(getLambdat.i(x.us),
      (row(matrix(0, nc, nc)) - 1)[upper.tri(matrix(0, nc, nc), diag = TRUE)])
    
    ## Testing getLambda in a different manner
    mat <- matrix(0, nc, nc)
    mat[upper.tri(mat, diag = TRUE)] <- getTheta(x.us)[getThetaIndex(x.us)]
    expect_equal(getLambda(x.us), t(mat))
   
    ## This is quite similar to what is in setMethod("getVC", ...) but using
    ## getPar() to make sure everything is working as planned
    ## TODO: think of a more creative method?
    if (nc <= 1L) {
      vcomp_test <- getPar(x.us)
      ccomp_test <- double(0L)
    } else {
      ii <- seq.int(from = 1L, by = nc + 1L, length.out = nc)
      i0 <- sequence.default(from = seq.int(from = 2L, by = nc + 1L, length.out = nc - 1L),
                             by = 1L,
                             nvec = (nc - 1L):1L)
      i1 <- sequence.default(from = ii,
                             by = 1L,
                             nvec = nc:1L)
      L <- matrix(0, nc, nc)
      L[i1] <- getPar(x.us)
      S <- tcrossprod(L)
      vcomp_test <- sqrt(S[ii])
      ccomp_test <- (S/vcomp_test/rep(vcomp_test, each = nc))[i0]
    }
    expect_equal(getVC(x.us)$vcomp, vcomp_test)
    expect_equal(getVC(x.us)$ccomp, ccomp_test)
  }
})

test_that("unit tests for diagonal covariances", {
  
  for (nc in c(0:4, 16L, 64L, 256L)){
    for(hom_test in c(TRUE, FALSE)){
      x.di <- new("Covariance.diag", nc = nc, hom = hom_test, simulate = TRUE)
      
      ## Basic sanity tests
      expect_s4_class(x.di, "Covariance.diag")
      expect_true(validObject(x.di))
      
      ## par should be the same as theta
      expect_equal(getPar(x.di), getTheta(x.di))
      
      ## checking for specific par/theta lengths
      if(hom_test){
        expect_equal(getParLength(x.di), if (nc > 0L) 1 else 0)
      } else {
        expect_equal(getParLength(x.di), if (nc > 0L) nc else 0)
      }
      
      expect_equal(getParLength(x.di), getThetaLength(x.di))
      expect_equal(length(getThetaIndex(x.di)), 
                   getThetaIndexLength(x.di))
      
      expect_equal(getLambdat.dp(x.di), rep(1L, nc))
      
      ## Testing getUpper() and getLower() in a different manner
      if(hom_test){
        expect_equal(getUpper(x.di), if(nc > 0L) Inf else numeric(0))
        expect_equal(getLower(x.di), if(nc > 0L) 0 else numeric(0))
      } else {
        expect_equal(getUpper(x.di), if(nc > 0L) rep(Inf, nc) else numeric(0))
        expect_equal(getLower(x.di), if(nc > 0L) rep(0, nc) else numeric(0))
      }
      
      ## Testing getLambdat.i in a different manner
      expect_equal(getLambdat.i(x.di), 
                   if(nc > 0L) 0:(nc-1) else integer(0))
      
      ## Testing getLambda in a different manner
      expect_equal(getLambda(x.di),
                   diag(nc) * getPar(x.di))
      
      ## Testing getVC
      expect_equal(getVC(x.di)$ccomp, numeric(0))
      expect_equal(getVC(x.di)$vcomp, getPar(x.di))
    }
  }
})

test_that("unit tests for compound symmetry covariances", {
  
  for (nc in c(0:4, 16L, 64L, 256L)){
    for(hom_test in c(TRUE, FALSE)){
      x.cs <- new("Covariance.cs", nc = nc, hom = hom_test, simulate = TRUE)
      
      ## Basic sanity tests
      expect_s4_class(x.cs, "Covariance.cs")
      expect_true(validObject(x.cs))
      
      ## par should be the same as vcomp and ccomp from getVC
      expect_equal(getPar(x.cs), c(getVC(x.cs)$vcomp, getVC(x.cs)$ccomp))
      
      ## testing for getTheta via the structure of getLambda
      ## (odd way of looking at things, but it nicely compares the relationship
      ## between these two functions.)
      lam_mat <- getLambda(x.cs)
      collected <- 1
      if(nc > 0L){
        if(hom_test){ ## for homogenous structures only
          collect <- numeric(nc*2-1)
          collected <- 1
          for(i in 1:nc){
            if(i != nc){
              collect[collected:(collected+1)] <- lam_mat[i:(i+1),i]
            } else {
              collect[collected:collected] <- lam_mat[i,i]
            }
            collected = collected + 2
          }
        } else { ## for non-homogenous structures
          grab <- nc
          collect <- numeric((nc * (nc+1)/2))
          for(i in 1:nc){
            if(i != nc){
              collect[collected:(collected + (nc-i))] <- lam_mat[i:nc,i]
            } else {
              collect[collected:collected] <- lam_mat[i,i]
            }
            collected = collected + (nc - i + 1)
          }
        }
      } else {
        collect = numeric(0)
      }
      expect_equal(getTheta(x.cs), collect)
      
      ## Testing getParLength
      if(hom_test){ ## for homogenous, should be no more than 2
        if(nc <= 2L){
          expect_equal(getParLength(x.cs), nc)
        } else {
          expect_equal(getParLength(x.cs), 2L)
        }
      } else { ## Comparing length to vcomp and ccomp 
        expect_equal(getParLength(x.cs),
                     sum(lengths(getVC(x.cs)[c("vcomp", "ccomp")])))
      }
      
      expect_equal(getParLength(x.cs), length(getPar(x.cs)))
      expect_equal(getThetaLength(x.cs), length(getTheta(x.cs)))
      
      ## Testing getThetaLength
      expect_equal(length(getThetaIndex(x.cs)), 
                   getThetaIndexLength(x.cs))
      if(hom_test){
        expect_equal(getThetaLength(x.cs), if(nc > 0L) (nc*2 - 1) else 0)
      } else {
        expect_equal(getThetaLength(x.cs), (nc * (nc+1)/2))
      }
      
      expect_equal(getLambdat.dp(x.cs), if(nc > 0L) 1:nc else integer(0))
      
      ## Testing getUpper() and getLower()
      if(nc > 1L){
        if(hom_test) {
          expect_equal(getUpper(x.cs), c(Inf, 1))
          expect_equal(getLower(x.cs), c(0.0, -1/(nc-1)))
        } else {
          expect_equal(getUpper(x.cs), c(rep(Inf, nc), 1))
          expect_equal(getLower(x.cs), c(rep(0.0, nc), -1/(nc-1)))
        }
      } else if (nc == 1L) {
        expect_equal(getUpper(x.cs), Inf)
        expect_equal(getLower(x.cs), 0)
      } else {
        expect_equal(getUpper(x.cs), numeric(0))
        expect_equal(getLower(x.cs), numeric(0))
      }
      
      ## Testing getLambdat.i in a different manner
      expect_equal(getLambdat.i(x.cs), 
                   (row(matrix(0, nc, nc)) - 1)[upper.tri(matrix(0, nc, nc), diag = TRUE)])
    }
  }
})

test_that("unit tests for autoregressive covariances", {
  
  for (nc in c(0:4, 16L, 64L, 256L)){
    for(hom_test in c(TRUE, FALSE)){
      x.ar1 <- new("Covariance.ar1", nc = nc, hom = hom_test, simulate = TRUE)
      
      ## Basic sanity tests
      expect_s4_class(x.ar1, "Covariance.ar1")
      expect_true(validObject(x.ar1))
      
      getVC_t <- getVC(x.ar1); vcomp_t <- getVC_t$vcomp; ccomp_t <- getVC_t$ccomp
      ## par should be the same as vcomp and ccomp from getVC
      expect_equal(getPar(x.ar1), c(vcomp_t, ccomp_t))
      
      ## getTheta test by relating to getVC
      v1 <- ccomp_t^(0L:(nc - 1L))
      v2 <- ccomp_t^(0L:(nc - 2L)) * sqrt(1 - ccomp_t^2)
      if(hom_test){
        if(nc > 1L){
          theta_test <- vcomp_t * c(v1, v2)
        } else {
          theta_test <- vcomp_t
        }
      } else {
        if(nc > 1L){
          theta_test <- vcomp_t[sequence.default(from = 1L:nc, nvec = nc:1L)] *
            c(v1, v2[sequence.default(from = 1L, nvec = (nc - 1L):1L)])
        } else {
          theta_test <- vcomp_t
        }
      }
      expect_equal(getTheta(x.ar1), theta_test)
      
      ## Testing for theta and par lengths
      expect_equal(getThetaLength(x.ar1), length(getTheta(x.ar1)))
      expect_equal(getParLength(x.ar1), length(getPar(x.ar1)))
      expect_equal(getParLength(x.ar1),
                   sum(lengths(getVC(x.ar1)[c("vcomp", "ccomp")])))
      
      lam_test <- matrix(0, nc, nc)
      if(nc > 0L){
        lam_test[lower.tri(lam_test, diag = TRUE)] <-
          if(hom_test){
            theta_test[sequence.default(from = rep(c(1L, nc + 1L), 
                                                   c(1L, nc - 1L)), nvec = nc:1L)]
          } else {
            theta_test
          }
      }
      expect_equal(getLambda(x.ar1), lam_test)
      
      expect_equal(getLambdat.dp(x.ar1), if(nc > 0L) 1:nc else integer(0))
      
      ## Testing getUpper and getLower
      if(nc > 1L){
        if(hom_test) {
          expect_equal(getUpper(x.ar1), c(Inf, 1))
          expect_equal(getLower(x.ar1), c(0, -1))
        } else {
          expect_equal(getUpper(x.ar1), c(rep(Inf, nc), 1))
          expect_equal(getLower(x.ar1), c(rep(0.0, nc), -1))
        }
      } else if (nc == 1L) {
        expect_equal(getUpper(x.ar1), Inf)
        expect_equal(getLower(x.ar1), 0)
      } else {
        expect_equal(getUpper(x.ar1), numeric(0))
        expect_equal(getLower(x.ar1), numeric(0))
      }
      
      ## Testing getLambdat.i in a different manner
      expect_equal(getLambdat.i(x.ar1), 
                   (row(matrix(0, nc, nc)) - 1)[upper.tri(matrix(0, nc, nc), diag = TRUE)])
    }
  }
})

other_mod <- readRDS("../../inst/testdata/test-covariance_structures_data.rds")

## Getting all equal as a number (in the all.equal examples documentation;
## don't know why they didn't make an argument instead!?)
## TODO: May want to move this to utilities...?
all.eqNum <- function(...) as.numeric(sub(".*:", '', all.equal(...)))

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
fm1.us <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy, REML = FALSE)
fm1.cs <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)
fm1.diag <- lmer(Reaction ~ Days + diag(Days | Subject), 
                 sleepstudy, REML = FALSE)

fm1.REML <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm1.us.REML <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
fm1.cs.REML <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
fm1.diag.REML <- lmer(Reaction ~ Days + diag(Days | Subject), 
                      sleepstudy)

## TODO: need to make sure we also have the ar1 versions once it's working
## properly...

test_that("integration tests for coef and fixef", {
  ## Ensuring fm1 gives the same result as fm1.us
  expect_equal(coef(fm1), coef(fm1.us))
  expect_equal(coef(fm1.REML), coef(fm1.us.REML))
  expect_equal(fixef(fm1), fixef(fm1.us))
  expect_equal(fixef(fm1.REML), fixef(fm1.us.REML))
  
  ## One of the expected summaries
  tmpf <- function(x) capture.output(print(summary(x),digits=1))
  tfun <- function(cc) {
    w <- grep("Fixed effects:", cc)
    cc[w:length(cc)]
  }
  expected_summary <- c("Fixed effects:",                         
                        "            Estimate Std. Error t value",
                        "(Intercept)      251          7      38",
                        "Days              10          2       7",
                        "",                                       
                        "Correlation of Fixed Effects:",          
                        "     (Intr)",                            
                        "Days -0.138")  
  expect_equal(tfun(tmpf(fm1)), expected_summary)
  
})

test_that("integration tests for sigma", {
  ## Ensuring fm1 gives the same result as fm1.us
  expect_equal(sigma(fm1), sigma(fm1.us))
  expect_equal(sigma(fm1.REML), sigma(fm1.us.REML))
  
  ## Ensuring computing sigma is consistent for lme4 
  ## against other packages
  expect_equal(all.eqNum(sigma(fm1), other_mod$fm1.glmmTMB_sigma),
               0, tol = 5e-5)
  expect_equal(all.eqNum(sigma(fm1.cs), other_mod$fm1.glmmTMB.cs_sigma),
               0, tol = 5e-5)
  expect_equal(all.eqNum(sigma(fm1.diag), other_mod$fm1.glmmTMB.diag_sigma),
               0, tol = 5e-5)
  expect_equal(all.eqNum(sigma(fm1), other_mod$fm1.nlme_sigma),
               0, tol = 5e-5)
  expect_equal(all.eqNum(sigma(fm1.cs), other_mod$fm1.nlme.cs_sigma),
               0, tol = 5e-5)
  expect_equal(all.eqNum(sigma(fm1.REML), other_mod$fm1.nlme.REML_sigma),
               0, tol = 5e-5)
  expect_equal(all.eqNum(sigma(fm1.cs.REML), other_mod$fm1.nlme.cs.REML_sigma),
               0, tol = 5e-5)
})

test_that("integration tests for vcov", {
  
  ## Ensuring fm1 gives the same result as fm1.us
  expect_equal(vcov(fm1), vcov(fm1.us))
  expect_equal(vcov(fm1.REML), vcov(fm1.us.REML))
  
  ## Ensuring variance-covariance matrix are consistent for lme4 
  ## against other packages
  expect_equal(all.eqNum(vcov(fm1), other_mod$fm1.glmmTMB_vcov,
                         check.class = FALSE,
                         check.attributes = FALSE),
               0, tol = 5e-5)
  
  expect_equal(all.eqNum(vcov(fm1.cs), 
                         other_mod$fm1.glmmTMB.cs_vcov,
                         check.class = FALSE,
                         check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(vcov(fm1.diag), 
                         other_mod$fm1.glmmTMB.diag_vcov,
                         check.class = FALSE,
                         check.attributes = FALSE),
               0, tol = 5e-5)
  
  expect_equal(all.eqNum(vcov(fm1), 
                         other_mod$fm1.nlme_vcov,
                         check.class = FALSE,
                         check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(vcov(fm1.cs), 
                         other_mod$fm1.nlme.cs_vcov,
                         check.class = FALSE,
                         check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(vcov(fm1.REML), 
                         other_mod$fm1.nlme.REML_vcov,
                         check.class = FALSE,
                         check.attributes = FALSE),
               0, tol = 5e-5)
  ## TODO: ugh, fix another issue...
  #expect_equal(all.eqNum(as.numeric(vcov(fm1.cs.REML)), 
  #                       as.numeric(other_mod$fm1.nlme.cs.REML_vcov)),
  #             0, tol = 5e-5)
  
})

test_that("integration tests for VarCorr", {
  
  ## Ensuring fm1 gives the same result as fm1.us
  expect_equal(VarCorr(fm1), VarCorr(fm1.us))
  expect_equal(VarCorr(fm1.REML), VarCorr(fm1.us.REML))

  ## Ensuring variance components are consistent for lme4 
  ## against other packages
  x1 <- c(as.matrix(VarCorr(fm1)[[1]]))
  x1.REML <- c(as.matrix(VarCorr(fm1.REML)[[1]]))
  y1 <- c(other_mod$fm1.glmmTMB_var)
  
  x2 <- c(as.matrix(VarCorr(fm1.us)[[1]]))
  y2 <- c(other_mod$fm1.glmmTMB.us_var)
  z2 <- c(other_mod$fm1.nlme_var)
  z2.REML <- c(other_mod$fm1.nlme.REML_var)
  
  x3 <- c(as.matrix(VarCorr(fm1.cs)[[1]]))
  x3.REML <- c(as.matrix(VarCorr(fm1.cs.REML)[[1]]))
  y3 <- c(other_mod$fm1.glmmTMB.cs_var)
  z3 <- c(other_mod$fm1.nlme.cs_var)
  z3.REML <- c(other_mod$fm1.nlme.cs.REML_var)
  
  x4 <- c(as.matrix(VarCorr(fm1.diag)[[1]]))
  y4 <- c(other_mod$fm1.glmmTMB.diag_var)
  
  ## Testing unstructured
  expect_equal(all.eqNum(x1, y1), 0, tol = 5e-5)
  expect_equal(all.eqNum(x2, y2), 0, tol = 5e-5)
  expect_equal(all.eqNum(x1.REML, z2.REML), 0, tol = 5e-5)
  ## TODO: below fails!
  #expect_equal(all.eqNum(x2, z2), 0, tol = 5e-5)
  ## Testing compound symmetry
  expect_equal(all.eqNum(x3, y3), 0, tol = 5e-5)
  expect_equal(all.eqNum(x3, z3), 0, tol = 5e-5)
  
  # TODO: WHY IS BELOW NOT MATCHING??? 
  #expect_equal(all.eqNum(x3.REML, z3.REML), 0, tol = 5e-5)
  ## Testing diag
  expect_equal(all.eqNum(x4, y4), 0, tol = 5e-5)
})

test_that("integration tests for ranef", {
  
  ## Ensuring fm1 gives the same result as fm1.us
  expect_equal(ranef(fm1), ranef(fm1.us))
  expect_equal(ranef(fm1.REML), ranef(fm1.us.REML))
  
  ## Ensuring extracting random modes of the random effects 
  ## are consistent for lme4 against other packages
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.glmmTMB_ranef),
                         as.matrix(ranef(fm1)$Subject)),
               0, tol = 5e-5)
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.glmmTMB.cs_ranef),
                         as.matrix(ranef(fm1.cs)$Subject)),
               0, tol = 5e-5)
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.glmmTMB.diag_ranef),
                         as.matrix(ranef(fm1.diag)$Subject)),
               0, tol = 5e-5)
  
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.nlme_ranef), 
                         as.matrix(ranef(fm1)$Subject)),
               0, tol = 5e-5)
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.nlme.REML_ranef), 
                         as.matrix(ranef(fm1.REML)$Subject)),
               0, tol = 5e-5)
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.nlme.cs_ranef), 
                         as.matrix(ranef(fm1.cs)$Subject)),
               0, tol = 5e-5)
  expect_equal(all.eqNum(as.matrix(other_mod$fm1.nlme.cs.REML_ranef), 
                         as.matrix(ranef(fm1.cs.REML)$Subject)),
               0, tol = 5e-5)
  
})

test_that("integration tests for predict", {
  
  ## Ensuring fm1 gives the same result as fm1.us
  expect_equal(predict(fm1), predict(fm1.us))
  expect_equal(predict(fm1.REML), predict(fm1.us.REML))
  
  ## Ensuring extracting predictions are consistent for 
  ## lme4 against other packages
  
  expect_equal(all.eqNum(other_mod$fm1.glmmTMB_predict,
                         predict(fm1), check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(other_mod$fm1.glmmTMB.cs_predict,
                         predict(fm1.cs), check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(other_mod$fm1.glmmTMB.diag_predict,
                         predict(fm1.diag), check.attributes = FALSE),
               0, tol = 5e-5)
  
  expect_equal(all.eqNum(other_mod$fm1.nlme_predict,
                         predict(fm1), check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(other_mod$fm1.nlme.REML_predict,
                         predict(fm1.REML), check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(other_mod$fm1.nlme.cs_predict,
                         predict(fm1.cs), check.attributes = FALSE),
               0, tol = 5e-5)
  expect_equal(all.eqNum(other_mod$fm1.nlme.cs.REML_predict,
                         predict(fm1.cs.REML), check.attributes = FALSE),
               0, tol = 5e-5)
})
