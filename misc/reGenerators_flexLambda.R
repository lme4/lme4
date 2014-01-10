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
					warning(
						paste("Using the experimental stuff when you could just use the stable specification?",
							  "What are you, some jerk?"))
				}
				
				cnms <- Ztl$cnms
				rm(Ztl)
				
				#<nc> variances +  1 corr
				ntheta <- nc + 1
				nlambda <- nl*((nc+1)*nc/2)
				
				# upper/lower limits:
				upper <- c(rep(Inf, nc), .99)
				lower <- c(rep(0, nc), -.99)
				
				#diagonal entries of the cholesky are sqrt(<this expression>)*sd
				diagfactors <- function(x, col=nc) -(((col-1)*x^2- (col-2)*x - 1) /((col-2)*x + 1))
				if(nc>2){
					#rational function is wild, so find small interval with root in it first
					grid <- seq(0, -1, l=500)
					if(any(diagfactors(grid)<0)){
						lower[nc+1] <- .95*uniroot(diagfactors,
                                                                           lower=grid[min(which(diagfactors(grid)<0))],
                                                                           upper=grid[min(which(diagfactors(grid)<0))-1])$root
					} 
				} 	
				
				if(!is.null(init)){
					stopifnot(length(init)!=nc+1)
					stopifnot(all(init[1:nc]>0))
					stopifnot((init[nc+1] > lower[nc+1]) & (init[nc+1] < upper[nc+1]))
					theta <- init 
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
					function(theta){
						#get rho
						rho <-  theta[nc+1]
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
				upper <- c(Inf, .99)
				lower <- c(0, -.99)
				
				#diagonal entries of the cholesky are sqrt(<this expression>)*sd
				diagfactors <- function(x, col=nc) -(((col-1)*x^2- (col-2)*x - 1) /((col-2)*x + 1))
				
				if(!is.null(init)){
					stopifnot(length(init)!=2)
					stopifnot(init[1]>0)
					stopifnot((init[2] > lower[2]) & (init[2] < upper[2]))
					theta <- init
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
					function(theta){
						#get rho
						rho <-  theta[2]
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


#' @title Autocorrelated random intercepts 
#' 
#' For a \code{formula=~(<time>|<id>)}, specifies random effects with an AR(1)-correlation
#' structure in (discrete, equidistant!) \code{<time>} for each level of \code{<id>}, i.e. for 2 observations
#' \eqn{y_ij} and \eqn{y_ij'} that are \eqn{d} units of \code{time} apart, the covariance is
#' \eqn{\sigma^2\rho^d}. 
#' In this first stab at this, observations within each group have to be ordered and 1 timeunit apart, because the math for
#' the entries in Lambdat gets too messy otherwise... 
#' The default, \code{formula=~(.|1)}, yields one auto-correlated intercept per observation under the assumptions
#' that the data are ordered and equidistant in time. Note that the covariance will be dense and slow as hell
#' to evaluate.
#' @param formula a one sided formula specifying a random effect in
#'     \code{lme4} notation, i.e. \code{~(<time> | <id>)}. \code{~(. | <id>)} is a valid
#'     specification that simply uses the rownumbers as time index. \code{<time>} cannot be a factor.
#' @param init (optional) initial values for the standard deviations and the correlation.
#' 		  If not supplied, sd's will be 1 and the inital correlation will be .1.	    
#' @param het NOT YET IMPLEMENTED
#' @param max.lag NOT YET IMPLEMENTED
#' @return a function creating a return object like \code{mkReTrm} 
ar1d <- function(formula=~(.|1), order=1, init=c(1, .2), het=NULL, max.lag=NULL){
	stopifnot(all(order==1))
	
	mkReTrmAR <- local({
		bar <- formula[[2]][[2]]
		addRowIndex <- deparse(bar[[2]])=="."
		
		function(fr){
	        
			ff <- getGrouping(bar, fr)
			nl <- length(levels(ff))
			if(addRowIndex){
				##FIXME: all this will surely break for predict/simulate...
				## still a nice default behaviour to have.
				fr[".rows."] <- 1:nrow(fr)
				Zt <-  Matrix(0, nrow(fr), nrow(fr))
				diag(Zt) <- 1
				Zt <- as(Zt, "dgCMatrix")
				nc <- NA ## makes no sense as can have differnt number of obs
				## per level...
				cnms <- ".rows."
				bar[[2]] <- as.symbol(".rows.")
			} else {
				#convert time variable into factor to get
				#one effect per timepoint per subject
				Zfr <- fr
				Zfr[[deparse(bar[[2]])]] <- as.factor(fr[[deparse(bar[[2]])]])
				#initialize transposed design
				Zbar <- bar
				Zbar[[2]] <- substitute( 0 + lhs, list(lhs=bar[[2]])  )
				Ztl <- mkZt(ff, Zbar, Zfr)
				Zt <- Ztl$Zt
				nc <- Ztl$nc
				cnms <- Ztl$cnms
				rm(Ztl)
			}
			#1 variance + 1 corr
			ntheta <- 2
			
			upper <- c(Inf, .99)
			lower <- c(0, -.99)
			
			#initalize Lambdatx
			arChol <- function(r, d, firstrow){
				firstrow*r^d + (!firstrow)*(sqrt(1-r^2))*r^d
			}
			timepoints <- with(fr, split(eval(bar[[2]]), ff))
			arinfo <- lapply(timepoints, function(t){
				#FIXME: we need only unique(t), but that
				#means the check for ordered, equidist is sloppy!
				t <- unique(t)
				if(any(diff(t)!=1)){
					stop("Timepoints either unsorted or not equidistant.")	
				}
				d <- abs(outer(t, t, "-"))
				firstrow <- as.vector((row(d)==1)[upper.tri(d, diag=TRUE)])
				d <- as.vector(d[upper.tri(d, diag=TRUE)])
				return(list(d=d, firstrow=firstrow, dim=length(t)))
			}) 
			Lambdablocks <- lapply(arinfo, function(info){
				fill <- arChol(init[2], info$d, info$firstrow)
				Lambdablock <- Matrix(0, info$dim, info$dim)
				Lambdablock[upper.tri(Lambdablock, diag=TRUE)] <- fill
				drop0(Lambdablock)
			})
			Lambdat <- init[1] * do.call(bdiag, unname(Lambdablocks))
			nlambda <- length(Lambdat@x)
			
			updateLambdatx <- local({
				# theta[1:nc]: sd's; theta[nc+1]: transformed corr
				arChol <- arChol
				d <- unlist(sapply(arinfo, "[[", "d"))
				firstrow <- unlist(sapply(arinfo, "[[", "firstrow"))
				function(theta){
					theta[1] * arChol(theta[2], d, firstrow)
				}
			})
			
			list(ff = ff, Zt = Zt, nl = nl, cnms = cnms,
				 nb = nrow(Zt), #how many ranefs
				 ntheta = ntheta, # how many var-cov. params
				 nc = nc, #how many ranefs per level
				 nlambda = nlambda, #how many non-zeroes in Lambdat 
				 Lambdat=Lambdat,
				 theta = init,
				 Lind = rep(NA, nlambda),
				 updateLambdatx = updateLambdatx,
				 upper = upper,
				 lower = lower, 
				 special = TRUE)
		}
		
	})
}

#' @title Random intercepts as realisations of a Gaussian random field with a given, 
#' fixed correlation structure. 
#' 
# grf <- function(formula=~(.|1), S){
# 	C <- chol(S)
# 	
# 	mkReTrmGrf <- local({
# 		C <- C
# 		bar <- formula[[2]][[2]]
# 		
# 		function(fr){
# 			
# 			ff <- getGrouping(bar, fr)
# 			nl <- length(levels(ff))
# 			stopifnot(is.factor(fr[[deparse(bar[[2]])]]))
# 			stopifnot(nlevels(fr[[deparse(bar[[2]])]]) == ncol(C))
# 			stopifnot(all(levels(fr[[deparse(bar[[2]])]]) == colnames(C)))
# 			
# 		    bar[[2]] <- substitute( 0 + lhs, list(lhs=bar[[2]])  )
# 			Ztl <- mkZt(ff, bar, fr)
# 			Zt <- Ztl$Zt
# 			nc <- Ztl$nc
# 			cnms <- Ztl$cnms
# 			rm(Ztl)
# 			
# 			#1 variance
# 			ntheta <- 1
# 			
# 			upper <- c(Inf)
# 			lower <- c(0)
# 			
# 			#initalize Lambdatx
# 			arChol <- function(r, d, firstrow){
# 				firstrow*r^d + (!firstrow)*(sqrt(1-r^2))*r^d
# 			}
# 			timepoints <- with(fr, split(eval(bar[[2]]), ff))
# 			arinfo <- lapply(timepoints, function(t){
# 				#FIXME: we need only unique(t), but that
# 				#means the check for ordered, equidist is sloppy!
# 				t <- unique(t)
# 				if(any(diff(t)!=1)){
# 					stop("Timepoints either unsorted or not equidistant.")	
# 				}
# 				d <- abs(outer(t, t, "-"))
# 				firstrow <- as.vector((row(d)==1)[upper.tri(d, diag=TRUE)])
# 				d <- as.vector(d[upper.tri(d, diag=TRUE)])
# 				return(list(d=d, firstrow=firstrow, dim=length(t)))
# 			}) 
# 			Lambdablocks <- lapply(arinfo, function(info){
# 				fill <- arChol(init[2], info$d, info$firstrow)
# 				Lambdablock <- Matrix(0, info$dim, info$dim)
# 				Lambdablock[upper.tri(Lambdablock, diag=TRUE)] <- fill
# 				drop0(Lambdablock)
# 			})
# 			Lambdat <- init[1] * do.call(bdiag, unname(Lambdablocks))
# 			nlambda <- length(Lambdat@x)
# 			
# 			updateLambdatx <- local({
# 				# theta[1:nc]: sd's; theta[nc+1]: transformed corr
# 				arChol <- arChol
# 				d <- unlist(sapply(arinfo, "[[", "d"))
# 				firstrow <- unlist(sapply(arinfo, "[[", "firstrow"))
# 				function(theta){
# 					theta[1] * arChol(theta[2], d, firstrow)
# 				}
# 			})
# 			
# 			list(ff = ff, Zt = Zt, nl = nl, cnms = cnms,
# 				 nb = nrow(Zt), #how many ranefs
# 				 ntheta = ntheta, # how many var-cov. params
# 				 nc = nc, #how many ranefs per level
# 				 nlambda = nlambda, #how many non-zeroes in Lambdat 
# 				 Lambdat=Lambdat,
# 				 theta = init,
# 				 Lind = rep(NA, nlambda),
# 				 updateLambdatx = updateLambdatx,
# 				 upper = upper,
# 				 lower = lower, 
# 				 special = TRUE)
# 		}
# 		
# 	})
# }

