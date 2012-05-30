## convert matrix list to concatenated vector of lower triangles
mlist2vec <- function(L) {
    n <- sapply(L,nrow)
    ## allow for EITHER upper- or lower-triangular input
    ff <- function(x) {
        if (all(x[upper.tri(x)]==0)) t(x[lower.tri(x,diag=TRUE)])
        else x[upper.tri(x,diag=TRUE)]
    }
    r <- unlist(lapply(L,ff))
    attr(r,"clen") <- n
    r
}

get_clen <- function(v,n=NULL) {
    if (is.null(n)) {
        if (is.null(n <- attr(v,"clen"))) {
            ## single component
            n <- (sqrt(8*length(v)+1)-1)/2
        }
    }
    n
}

## convert concatenated vector to matrix list
## (lower triangle or symmetric)
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

## convert 'sdcor' format -- diagonal = std dev, off-diag=cor
##  to and from variance-covariance matrix
sdcor2cov  <- function(m) {
    sd <- diag(m)
    diag(m) <- 1
    m * outer(sd,sd)
}

cov2sdcor  <- function(m) {
    v <- diag(m)
    m1 <- cov2cor(m)
    diag(m1) <- sqrt(v)
    m1
}

dmult <- function(m,s) {
    diag(m) <- diag(m)*s
    m
}

safe_chol <- function(m) {
    if (all(m==0)) return(m)
    if (nrow(m)==1) return(sqrt(m))
    if (all(dmult(m,0)==0)) {  ## diagonal
        return(diag(sqrt(diag(m))))
    }
    ## attempt regular Chol. decomp
    if (!inherits(try(cc <- chol(m),silent=TRUE),"try-error"))
        return(cc)
    ## ... pivot if necessary ...
    cc <- suppressWarnings(chol(m,pivot=TRUE))
    oo <- order(attr(cc,"pivot"))
    cc[,oo]
    ## FIXME: pivot is here to deal with semidefinite cases,
    ## but results might be returned in a strange format: TEST
}

Vv_to_Cv <- function(v,n=NULL,s=1) {
    if (!missing(s)) {
        v <- v[-length(v)]
    }
    r <- mlist2vec(lapply(vec2mlist(v,n,symm=TRUE),
                          function(m) t(safe_chol(m/s^2))))
    attr(r,"clen") <- get_clen(v,n)
    r
}

Sv_to_Cv <- function(v,n=NULL,s=1) {
    if (!missing(s)) {
        v <- v[-length(v)]
    }
    r <-  mlist2vec(lapply(vec2mlist(v,n,symm=TRUE),
                           function(m) t(safe_chol(sdcor2cov(m)/s^2))))
    attr(r,"clen") <- get_clen(v,n)
    r
}

Cv_to_Vv <- function(v,n=NULL,s=1) {
    r <- mlist2vec(lapply(vec2mlist(v,n,symm=FALSE),
                          function(m) tcrossprod(m)*s^2))
    if (!missing(s)) r <- c(r,s^2)
    attr(r,"clen") <- get_clen(v,n)
    r
}

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

