#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1,v2){-log(pi)-log(v2)-log(1+((x-v1)/v2)^2)}
cauchy_logfdd=Deriv(f,c("v1","v2"),nderiv=2)
cauchy_logfddd=Deriv(f,c("v1","v2"),nderiv=3)
#
sink("42c_cauchy_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("cauchy_logfdd=")
print.function(cauchy_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("cauchy_logfddd=")
print.function(cauchy_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"cauchy_ldda=function(x,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(cauchy_logfdd)
	temp1=vf(x,v1,v2)
	temp2=apply(temp1,1,sum)/nx
	for (i in 1:2){
		for (j in 1:2){
			ldd[i,j]=temp2[(i-1)*2+j]
		}
	}
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"cauchy_lddda=function(x,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(cauchy_logfddd)
	temp1=vf(x,v1,v2)
	temp2=apply(temp1,1,sum)/nx
	for (i in 1:2){
		for (j in 1:2){
			for (k in 1:2){
				lddd[i,j,k]=temp2[(i-1)*4+(j-1)*2+k]
			}
		}
	}
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
