#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1,v2,v3){-log(v2)-(1+1/v3)*log(1+v3*(x-v1)/v2)-(1+v3*(x-v1)/v2)^(-1/v3)}
gev_logfdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
gev_logfddd=Deriv(f,c("v1","v2","v3"),nderiv=3)
#
sink("110c_gev_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_logfdd=")
print.function(gev_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_logfddd=")
print.function(gev_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_ldda=function(x,v1,v2,v3){
	nx=length(x)
	ldd=matrix(0,3,3)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfdd)
	temp1=vf(x,v1,v2,v3)
	temp2=apply(temp1,1,sum)/nx
	for (i in 1:3){
		for (j in 1:3){
			ldd[i,j]=temp2[(i-1)*3+j]
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
"gev_lddda=function(x,v1,v2,v3){
	nx=length(x)
	lddd=array(0,c(3,3,3))

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfddd)
	temp1=vf(x,v1,v2,v3)
	temp2=apply(temp1,1,sum)
	for (i in 1:3){
		for (j in 1:3){
			for (k in 1:3){
				lddd[i,j,k]=temp2[(i-1)*9+(j-1)*3+k]
			}
		}
	}
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
