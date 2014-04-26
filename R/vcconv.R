## These files are not currently exported; we are still trying to figure
## out the most appropriate user interface/naming convention/etc
## These functions will become more important, and need to be
## modified/augmented, if we start allowing for varying V-C structures
## (e.g. diagonal matrices, compound symmetry ...)

## In principle we might want to extract or input information:

## 1. as variance-covariance matrices
## 2. as Cholesky factors
## 3. as 'sdcorr' matrices (std dev on diagonal, correlations off diagonal)

## and we might want the structure to be:

## 1. a concatenated vector representing the lower triangles
##    (with an attribute carrying the information about group sizes)
## 2. a list of lower-triangle vectors
## 3. a list of matrices
## 4. a block-diagonal matrix

## If we are trying to convert to and from theta vectors, we also
## have to consider whether we are returning scaled Cholesky factors/
## var-cov matrices or unscaled ones.  For the code below I have
## chosen to allow the residual variance etc. to be appended as
## the last element of a variance-covariance vector.  (This last
## part is a little less generic than the rest of it.)

##' List of matrices to concatenated vector
##'
##' Convert list of matrices to concatenated vector of lower triangles
##' with an attribute that gives the dimension of each matrix in the
##' original list.  This attribute may be used to reconstruct the
##' matrices.
##'
##' @param L list of symmetric, upper-triangular, or lower-triangular
##' square matrices
##' @return A concatenation of the elements in one triangle of each
##' matrix. An attribute \code{"clen"} gives the dimension of each
##' matrix.
mlist2vec <- function(L) {
    if(is.atomic(L)) L <- list(L)
    n <- sapply(L,nrow)
    ## allow for EITHER upper- or lower-triangular input;
    ## in either case, read off in "lower-triangular" order
    ## (column-wise)
    ff <- function(x) {
        if (all(x[upper.tri(x)]==0)) t(x[lower.tri(x,diag=TRUE)])
        else t(x)[lower.tri(x,diag=TRUE)]
    }
    r <- unlist(lapply(L,ff))
    attr(r,"clen") <- n
    r
}

## Compute dimensions of a square matrix from the size
## of the lower triangle (length as a vector)
get_clen <- function(v,n=NULL) {
    if (is.null(n)) {
        if (is.null(n <- attr(v,"clen"))) {
            ## single component
            n <- (sqrt(8*length(v)+1)-1)/2
        }
    }
    n
}

##' Concatenated vector to list of matrices
##' 
##' Convert concatenated vector to list of matrices (lower triangle or
##' symmetric). These matrices could represent Cholesky factors,
##' covariance matrices, or correlation matrices (with standard
##' deviations on the diagonal).
##'
##' @param v concatenated vector
##' @param n FIXME: this has something to do with the dimension of
##' associated matrices.
##' @param symm Return symmetric matrix if \code{TRUE} or
##' lower-triangular if \code{FALSE}
##' @return List of matrices
vec2mlist <- function(v,n=NULL,symm=TRUE) {
    n <- get_clen(v,n)
    s <- split(v,rep.int(seq_along(n),n*(n+1)/2))
    m <- mapply(function(x,n0) {
        m0 <- diag(nrow=n0)
        m0[lower.tri(m0,diag=TRUE)] <- x
        if (symm) m0[upper.tri(m0)] <- t(m0)[upper.tri(m0)]
        m0
    },s,n,SIMPLIFY=FALSE)
    m
}

## Convert concatenated vector to list of ST matrices
vec2STlist <- function(v, n = NULL){
  ch <- vec2mlist(v, n, FALSE) # cholesky
  nch <- length(ch)
  sdiag <- function(x) { ## 'safe' diag()
      if (length(x)==1) matrix(x,1,1) else diag(x)
  }
  lapply(ch, function(L) {
    ST <- L%*%sdiag(1/sdiag(L))
    diag(ST) <- diag(L)
    ST
  })
}

##' Standard deviation-correlation matrix to covariance matrix
##' 
##' convert 'sdcor' format -- diagonal = std dev, off-diag=cor to and
##' from variance-covariance matrix
##'
##' @param m Standard deviation-correlation matrix
##' @return Covariance matrix
sdcor2cov  <- function(m) {
    sd <- diag(m)
    diag(m) <- 1
    m * outer(sd,sd)
}

##' Covariance matrix to standard deviation-correlation matrix
##'
##' convert cov to sdcor
##'
##' @param V Covariance matrix
##' @return Standard deviation-correlation matrix
cov2sdcor  <- function(V) {
    ## "own version" of cov2cor():  1. no warning for NA;   2. diagonal = sd(.)
    p <- (d <- dim(V))[1L]
    if (!is.numeric(V) || length(d) != 2L || p != d[2L]) 
        stop("'V' is not a square numeric matrix")
    sd <- sqrt(diag(V))
    Is <- 1/sd
    r <- V
    r[] <- Is * V * rep(Is, each = p)
    diag(r) <- sd
    r
}

## dmult <- function(m,s) {
##     diag(m) <- diag(m)*s
##     m
## }
## From Matrix package  isDiagonal(.) :

## for pre-3.0.0 compatibility:
if (!exists("rep_len")) rep_len <- function(x, length.out) {
    rep(x,length.out=length.out)
}
                                            
all0 <- function(x) !any(is.na(x)) && all(!x)
.isDiagonal.sq.matrix <- function(M, n = dim(M)[1L])
    all0(M[rep_len(c(FALSE, rep.int(TRUE,n)), n^2)])

## attempt to compute Cholesky, allow for positive semi-definite cases
##  (hackish)
safe_chol <- function(m) {
    if (all(m==0)) return(m)
    if (nrow(m)==1) return(sqrt(m))
    if (.isDiagonal.sq.matrix(m)) return(diag(sqrt(diag(m))))

    ## attempt regular Chol. decomp
    if (!is.null(cc <- tryCatch(chol(m), error=function(e) NULL)))
	return(cc)
    ## ... pivot if necessary ...
    cc <- suppressWarnings(chol(m,pivot=TRUE))
    oo <- order(attr(cc,"pivot"))
    cc[,oo]
    ## FIXME: pivot is here to deal with semidefinite cases,
    ## but results might be returned in a strange format: TEST
}

##' Variance-covariance to relative covariance factor
##' 
##' from var-cov to scaled Cholesky:
##'
##' @param v Vector of elements from the lower triangle of a
##' variance-covariance matrix.
##' @param n FIXME: see "@param n" above
##' @param s Scale parameter
##' @return Vector of elements from the lower triangle of a relative
##' covariance factor. 
Vv_to_Cv <- function(v,n=NULL,s=1) {
    if (!missing(s)) {
        v <- v[-length(v)]
    }
    r <- mlist2vec(lapply(vec2mlist(v,n,symm=TRUE),
                          function(m) t(safe_chol(m/s^2))))
    attr(r,"clen") <- get_clen(v,n)
    r
}

##' Standard-deviation-correlation to relative covariance factor
##' 
##' from sd-cor to scaled Cholesky
##'
##' @param v Vector of elements from the lower triangle of a
##' standard-deviation-correlation matrix.
##' @param n FIXME: see "@param n" above
##' @param s Scale parameter
##' @return Vector of elements from the lower triangle of a relative
##' covariance factor. 
Sv_to_Cv <- function(v,n=NULL,s=1) {
    if (!missing(s)) {
        v <- v[-length(v)]
    }
    r <-  mlist2vec(lapply(vec2mlist(v,n,symm=TRUE),
                           function(m) t(safe_chol(sdcor2cov(m)/s^2))))
    attr(r,"clen") <- get_clen(v,n)
    r
}

##' Relative covariance factor to variance-covariance
##' 
##' from unscaled Cholesky vector to (possibly scaled)
##' variance-covariance vector
##'
##' @param v Vector of elements from the lower triangle of a
##' relative covariance factor.
##' @param n FIXME: see "@param n" above
##' @param s Scale parameter
##' @return Vector of elements from the lower triangle of a
##' variance-covariance matrix.
Cv_to_Vv <- function(v,n=NULL,s=1) {
    r <- mlist2vec(lapply(vec2mlist(v,n,symm=FALSE),
                          function(m) tcrossprod(m)*s^2))
    if (!missing(s)) r <- c(r,s^2)
    attr(r,"clen") <- get_clen(v,n)
    r
}

##' Relative covariance factor to standard-deviation-correlation
##' 
##' from unscaled Chol to sd-cor vector
##'
##' @param v Vector of elements from the lower triangle of a
##' relative covariance factor.
##' @param n FIXME: see "@param n" above
##' @param s Scale parameter
##' @return Vector of elements from the lower triangle of a
##' standard-deviation-correlation matrix.
Cv_to_Sv <- function(v,n=NULL,s=1) {
    r <- mlist2vec(lapply(vec2mlist(v,n,symm=FALSE),
                          function(m) cov2sdcor(tcrossprod(m)*s^2)))
    if (!missing(s)) r <- c(r,s)
    attr(r,"clen") <- get_clen(v,n)
    r
}

if (FALSE) {
    cvec1 <- 1:6
    Cv_to_Vv(cvec1)
    Vv_to_Cv(Cv_to_Vv(0))
    Cv_to_Vv(cvec1,s=2)
    Sv_to_Cv(Cv_to_Sv(cvec1))
    Vv_to_Cv(Cv_to_Vv(cvec1))
    ## for length-1 matrices, Cv_to_Sv should be equivalent
    ##   to multiplying Cv by sigma and appending sigma ....
    clist2 <- list(matrix(1),matrix(2),matrix(3))
    cvec2 <- mlist2vec(clist2)
    all((cvec3 <- Cv_to_Sv(cvec2,s=2))==c(cvec2*2,2))
    all(Sv_to_Cv(cvec3,n=rep(1,3),s=2)==
        cvec3[-length(cvec3)]/cvec3[length(cvec3)])
}

