
test_that("unit test for unstructured covariances", {
  
  for(i in seq_len(20L)){
    nc <- sample(0:10, 1L)
    
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
  
  for(i in seq_len(20L)){
    nc <- sample(0:10, 1L)
    hom_test <- sample(c(T, F), size = 1)
    
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
})

test_that("unit tests for compound symmetry covariances", {
  
  for(i in seq_len(20L)){
    nc <- sample(0:10, 1L)
    hom_test <- sample(c(T, F), size = 1)
    
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
    # TODO: the test below currently fails
    #expect_equal(getThetaLength(x.cs), length(getTheta(x.cs)))
    
    ## Testing getThetaLength
    expect_equal(length(getThetaIndex(x.cs)), 
                 getThetaIndexLength(x.cs))
    if(hom_test){
      # TODO: the test below currently fails
      #expect_equal(getThetaLength(x.cs), if(nc > 0L) (nc*2 - 1) else 0)
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
})

test_that("unit tests for autoregressive covariances", {
  
  for(i in seq_len(20L)){
    nc <- sample(0:10, 1L)
    hom_test <- sample(c(T, F), size = 1)
    
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
})
