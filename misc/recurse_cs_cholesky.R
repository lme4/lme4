## recursive solution to Cholesky factor of compound symmetric matrix
sym_cschol <- function(n, rho) {
  L <- matrix(0, n, n)
  
  L[1, 1] <- 1
  
  for (i in 2:n) {
    ## Compute off-diagonal elements in row i
    for (j in 1:(i - 1)) {
      sum_prod <- sum(L[i, 1:(j - 1)] * L[j, 1:(j - 1)])
      L[i, j] <- (rho - sum_prod) / L[j, j]
    }
    
    ## Compute diagonal element in row i
    sum_squares <- sum(L[i, 1:(i - 1)]^2)
    L[i, i] <- sqrt(1 - sum_squares)
  }
  
  return(L)
}

num_cschol <- function(n, rho) {
  L <- matrix(rho, n, n)
  diag(L) <- 1
  t(chol(L))
}

all.equal(sym_cschol(5, 0.3), num_cschol(5, 0.3))
