#' @title Random effects with diagonal covariance
#' 
#' Specifies random effects without correlation, i.e., the covariance
#' for the \eqn{q}-dim. random effect in each grouping level is 
#' \eqn{\text{diag}(\vartheta^2_1, \dots, \vartheta^2_q))} if \code{iid==FALSE} or
#' \eqn{\vartheta^2 I_q} if \code{iid==TRUE}.
#' 
#' @param formula a one sided formula specifying a random effect in
#'     \code{lme4} notation, i.e. \code{~(<covariates> | <grouping>)}.
#' @param iid enforce identical variances for each component of the random effects.
#' @return a function creating a return object like \code{mkReTrm}    
d <- function(formula, iid=FALSE){
	mkReTrmDiagonal <- local({
		bar <- formula[[2]][[2]]
		iid <- iid
		function(fr){
			ff <- getGrouping(bar, fr)
			
			nl <- length(levels(ff))
			
			#initialize transposed design
			Ztl <- mkZt(ff, bar, fr)
			Zt <- Ztl$Zt
			nc <- Ztl$nc
			cnms <- Ztl$cnms
			rm(Ztl)
			
			#enforce identical variance for all effects?
			if(iid){
				ntheta <- 1 #how many var/cov-parameters
				Lind <- rep(1, times=nl*nc)
			} else {
				ntheta <- nc #how many var/cov-parameters
				## which theta goes where in Lambdat
				Lind <- rep(1:ntheta, times=nl)
			}
			
			#initialize variances with  1
			theta <- rep(1, ntheta)
			
			# initialize the upper triangular Cholesky factor for Cov(b)
			# we want standard dgC not ddI (?) but no direct coerce exists:
			Lambdat <- as(as(kronecker(Diagonal(nc), Diagonal(nl)), 
							 "dgTMatrix"), "dgCMatrix")
			
			
			# upper/lower limits:
			upper <- rep(Inf, ntheta)
			lower <- rep(0, ntheta)
			
			list(ff = ff, Zt = Zt, nl = nl, cnms = cnms,
				 nb = nl*nc, #how many ranefs
				 ntheta = ntheta, # how many var-cov. params
				 nc = nc, #how many ranefs per level
				 nlambda = nl*nc, #how many non-zeroes in Lambdat 
				 Lambdat=Lambdat,
				 theta = theta,
				 Lind = Lind,
				 updateLambdatx = local({
				 	Lind <- Lind
				 	function(theta) theta[Lind]
				 }),
				 upper = upper,
				 lower = lower, 
				 special = TRUE)
		}
	})
	mkReTrmDiagonal
}





ar <- function(formula=~(.rows|1), order=1, max.dist=NA, init=.2){
	stopifnot(all(order==1, band==NA))
	
	mkIndex <- if(deparse(bar[[2]])==".rows"){
		mk
	}
	
	mkReTrmAR <- local({
		bar <- formula[[2]][[2]]
		function(fr){
			ff <- as.factor(fr[[deparse(group[[2]])]])
			nl <- length(levels(ff))
			
			
			addRowIndex <- deparse(bar[[2]])==".rows"
			if(addRowIndex) fr[".rows"] <- 1:nrow(fr)
			
			
		}
	})
	
}