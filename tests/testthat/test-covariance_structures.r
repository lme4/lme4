
test_that("unit test for unstructured covariances", {
  
  for(i in seq_len(100L)){
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
