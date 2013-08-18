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
				ntheta <- 1 
				# which theta goes where in Lambdat
				Lind <- rep(1, times=nl*nc)
			} else {
				ntheta <- nc 
				# which theta goes where in Lambdat
				Lind <- rep(1:ntheta, times=nl)
			}
			
			#initialize variances with  1
			theta <- rep(1, ntheta)
			
			# initialize the upper triangular Cholesky factor for Cov(b)
			Lambdat <- as(do.call(bdiag, replicate(nl, list(Diagonal(nc)))),
						  "dgCMatrix")
			
			
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

#' @title (Variance-hetergoneous) compound symmetry random effects 
#' 
#' Specifies correlated random effects with constant correlation > 0, i.e., the covariance
#' between entries \eqn{b_{it}} and \eqn{b_{ir}} in the \eqn{q}-dim. random effect \eqn{b_i} 
#' for level \eqn{i} of the grouping factor is \eqn{\rho \sigma_r \sigma_t} (instead of
#' \eqn{\rho_{rt}\sigma_r\sigma_t} for the standard unstructured random effects). 
#' If \code{het=FALSE},  enforces variance homogeneity (\eqn{\sigma_r = \sigma_t \forall r,t}).
#' 
#' @param formula a one sided formula specifying a random effect in
#'     \code{lme4} notation, i.e. \code{~(<covariates> | <grouping>)}.
#' @param init (optional) initial values for the standard deviations and the correlation.
#' 		  If not supplied, sd's will be 1 and the inital correlation will be .1.	    
#' @param het allow variance heterogeneity? defaults to true.
#' @return a function creating a return object like \code{mkReTrm}    
cs <- function(formula, init=NULL, het=TRUE){
	if(het){
		mkReTrmCSHet <- local({
			bar <- formula[[2]][[2]]
			function(fr){
				ff <- getGrouping(bar, fr)
				
				nl <- length(levels(ff))
				
				#initialize transposed design
				Ztl <- mkZt(ff, bar, fr)
				Zt <- Ztl$Zt
				nc <- Ztl$nc
				
				if(nc <= 2){
					warning("Using the experimental stuff when you could just use the stable specification?
						What are you, som?")
				}
				
				cnms <- Ztl$cnms
				rm(Ztl)
				
				#<nc> variances +  1 corr
				ntheta <- nc + 1
				nlambda <- nl*((nc+1)*nc/2)
				
				# upper/lower limits:
				upper <- rep(Inf, nc+1)
				lower <- c(rep(0, nc), -Inf)
				
				#diagonal entries of the cholesky are sqrt(<this expression>)*sd
				diagfactors <- function(x, col=nc) -(((col-1)*x^2- (col-2)*x - 1) /((col-2)*x + 1))
				# (upper not implemented yet, so need transform on corr)
				# lower limit for corr depends on nc:
				# - ((nc-1)^2*rho-(nc-2)*rho - 1) / ((nc-2)*rho + 1) > 0
				cl <- c(-.5, .5) 
				if(nc>2){
					#rational function is wild, so find small interval with root in it first
					grid <- seq(0, cl[1], l=500)
					if(any(diagfactors(grid)<0)){
						cl[1] <- uniroot(diagfactors, lower=grid[min(which(diagfactors(grid)<0))],
										 upper=grid[min(which(diagfactors(grid)<0))-1])$root
					} 
				} 	
				
				# initialize variances (default sd=1, cor=<middle of valid range>)
				# (upper not implemented yet, so map reals into (cl[1], cl[2]) by a scaled logistic
				if(!is.null(init)){
					stopifnot(length(init)!=nc+1)
					stopifnot(all(init[1:nc]>0))
					stopifnot((init[nc+1] > cl[1]) & (init[nc+1] < cl[2]))
					theta <- init
					theta[nc+1] <- qlogis((init[nc+1]-cl[1])/(diff(cl))) 
				} else {
					theta <- c(rep(1, nc), .1)
				}
				
				
				# initialize the upper triangular Cholesky factor for Cov(b)
				Lambdat <- 1*bdiag(replicate(nl, 
											 list(upper.tri(diag(nc), diag=TRUE))))
				# somebody more Asperger than me can write up
				# the explicit analytic recursion for the entries 
				# in the chol-factor -- it's just one (nc x nc) Cholesky
				# per update, so brute-forcing it should be fine...
				updateLambdatx <- local({
					# theta[1:nc]: sd's; theta[nc+1]: transformed corr
					C <- diag(nc)
					nl <- nl
					cl <- cl
					function(theta){
						#get rho
						rho <-  diff(cl) * plogis(theta[nc+1]) + cl[1]
						C[upper.tri(C)] <- C[lower.tri(C)] <- rho
						#make covariance:
						sd <- theta[1:nc]
						S <- t(C*sd)*sd
						#drop zero rows/columns, take chol
						d0 <- which(sd < .Machine$double.eps)
						if(length(d0)){
							aux <- chol(S[-d0,-d0])
							cS <- diag(nc)
							diag(cS)[d0] <- 0
							cS[-d0, -d0] <- aux
						} else {
							cS <- chol(S)
						}	
						lambdatx1 <- as.vector(cS)[upper.tri(cS, diag=TRUE)]
						rep(lambdatx1, times=nl)
					}
				})
				Lambdat@x <- updateLambdatx(theta)	
				
				
				list(ff = ff, Zt = Zt, nl = nl, cnms = cnms,
					 nb = nl*nc, #how many ranefs
					 ntheta = ntheta, # how many var-cov. params
					 nc = nc, #how many ranefs per level
					 nlambda = nlambda, #how many non-zeroes in Lambdat 
					 Lambdat=Lambdat,
					 theta = theta,
					 Lind = rep(NA, nlambda),
					 updateLambdatx = updateLambdatx,
					 upper = upper,
					 lower = lower, 
					 special = TRUE)
			}
		})	
		return(mkReTrmCSHet)
	} else {
		mkReTrmCSHom <- local({
			bar <- formula[[2]][[2]]
			function(fr){
				ff <- getGrouping(bar, fr)
				
				nl <- length(levels(ff))
				
				#initialize transposed design
				Ztl <- mkZt(ff, bar, fr)
				Zt <- Ztl$Zt
				nc <- Ztl$nc
				
				cnms <- Ztl$cnms
				rm(Ztl)
				
				#1 variance + 1 corr
				ntheta <- 2
				nlambda <- nl*((nc+1)*nc/2)
				
				# upper/lower limits:
				upper <- rep(Inf, 2)
				lower <- c(0, -Inf)
				
				#diagonal entries of the cholesky are sqrt(<this expression>)*sd
				diagfactors <- function(x, col=nc) -(((col-1)*x^2- (col-2)*x - 1) /((col-2)*x + 1))
				# (upper not implemented yet, so need transform on corr)
				# lower limit for corr depends on nc:
				# - ((nc-1)^2*rho-(nc-2)*rho - 1) / ((nc-2)*rho + 1) > 0
				cl <- c(-.5, .5) 
				if(nc>2){
					#rational function is wild, so find small interval with root in it first
					grid <- seq(0, cl[1], l=500)
					if(any(diagfactors(grid)<0)){
						cl[1] <- uniroot(diagfactors, lower=grid[min(which(diagfactors(grid)<0))],
										 upper=grid[min(which(diagfactors(grid)<0))-1])$root
					} 
				} 	
				
				# initialize variances (default sd=1, cor=<middle of valid range>)
				# (upper not implemented yet, so map reals into (cl[1], cl[2]) by a scaled logistic
				if(!is.null(init)){
					stopifnot(length(init)!=2)
					stopifnot(init[1]>0)
					stopifnot((init[2] > cl[1]) & (init[2] < cl[2]))
					theta <- init
					theta[2] <- qlogis((init[2]-cl[1])/(diff(cl))) 
				} else {
					theta <- c(1, .1)
				}
				
				
				# initialize the upper triangular Cholesky factor for Cov(b)
				Lambdat <- 1*bdiag(replicate(nl, 
											 list(upper.tri(diag(nc), diag=TRUE))))
				# somebody more Asperger than me can write up
				# the explicit analytic recursion for the entries 
				# in the chol-factor -- it's just one (nc x nc) Cholesky
				# per update, so brute-forcing it should be fine...
				updateLambdatx <- local({
					# theta[1:nc]: sd's; theta[nc+1]: transformed corr
					C <- diag(nc)
					nl <- nl
					cl <- cl
					function(theta){
						#get rho
						rho <-  diff(cl) * plogis(theta[2]) + cl[1]
						C[upper.tri(C)] <- C[lower.tri(C)] <- rho
						#make covariance:
						sd <- theta[1]
						S <- t(C*sd)*sd
						#drop zero rows/columns, take chol
						d0 <- which(sd < .Machine$double.eps)
						if(length(d0)){
							aux <- chol(S[-d0,-d0])
							cS <- diag(nc)
							diag(cS)[d0] <- 0
							cS[-d0, -d0] <- aux
						} else {
							cS <- chol(S)
						}	
						lambdatx1 <- as.vector(cS)[upper.tri(cS, diag=TRUE)]
						rep(lambdatx1, times=nl)
					}
				})
				Lambdat@x <- updateLambdatx(theta)	
				
				
				list(ff = ff, Zt = Zt, nl = nl, cnms = cnms,
					 nb = nl*nc, #how many ranefs
					 ntheta = ntheta, # how many var-cov. params
					 nc = nc, #how many ranefs per level
					 nlambda = nlambda, #how many non-zeroes in Lambdat 
					 Lambdat=Lambdat,
					 theta = theta,
					 Lind = rep(NA, nlambda),
					 updateLambdatx = updateLambdatx,
					 upper = upper,
					 lower = lower, 
					 special = TRUE)
			}
			
		})
		return(mkReTrmCSHom)
	}
}
# 
# # fixed correlation
# fc <- function(formula, P, S=solve)
# 	
# 	fr <- expand.grid(t=1:5, subject=factor(1:10))
# formula <- ~(.|1)
# mkIndex <- if(deparse(bar[[2]])=="."){
# 	fr[["..rows"]] <- 1:nrow(fr)
# 	
# 	
# }
# 
# 
# 
# ar <- function(formula=~(.|1), order=1, max.dist=NA, init=.2){
# 	stopifnot(all(order==1, band==NA))
# 	
# 	mkIndex <- if(deparse(bar[[2]])=="."){
# 		mk
# 	}
# 	
# 	mkReTrmAR <- local({
# 		bar <- formula[[2]][[2]]
# 		function(fr){
# 			ff <- as.factor(fr[[deparse(group[[2]])]])
# 			nl <- length(levels(ff))
# 			
# 			
# 			addRowIndex <- deparse(bar[[2]])==".rows"
# 			if(addRowIndex) fr[".rows"] <- 1:nrow(fr)
# 			
# 			
# 		}
# 	})
# 	
# }