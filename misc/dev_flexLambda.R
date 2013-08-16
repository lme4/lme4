#setwd("C:/lme4")
#library(devtools)
#dev_mode() 
#load_all()
rm(list=ls())
#reload()


(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy))
(fm3 <- lmer(Reaction ~ Days + (0+ Days | Subject) + (1|Subject), sleepstudy))



n <- 100
data <- data.frame(y=rnorm(n), x1=rnorm(n), x2=rnorm(n),
				   id=gl(20, n/20), id2=sample(gl(10, n/10)))
reGenerators=NULL
frml <- y ~ x1 + (x1|id) + (x1|id2)

fr <- model.frame(formula=as.formula(subbars(frml), env=data),
				  data=data, drop.unused.levels=TRUE)
attr(fr,"formula") <- frml
bars <- findbars(frml[[3]])


ar <- function(formula=~(.rows|1), order=1, init=0){
	
	generator <- local({
		call <- match.call()	
		function(fr){
		    groupingfrml <- substitute(~ 0 + grp, 
		    						   list(grp=call$formula[[2]][[2]][[3]]))
		    blocks <- if(call$formula[[2]][[2]][[3]]==1)
		    						   
		}
	}) 
		
	
}


arGenerator <- function(fr){
	
}

reTrm <- mkReList(bars[[1]])
