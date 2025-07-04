#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t1,t2,v1,v2,v3,v4,v5){-v3-t2*v4-(1+1/v5)*log(1+v5*(x-v1-t1*v2)/exp(v3+t2*v4))-
		(1+v5*(x-v1-t*v2)/exp(v3+t2*v4))^(-1/v5)}
gev_p12_logfdd=Deriv(f,c("v1","v2","v3","v4","v5"),nderiv=2)
gev_p12_logfddd=Deriv(f,c("v1","v2","v3","v4","v5"),nderiv=3)
#
sink("151c_gev_p12_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_logfdd=")
print.function(gev_p12_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_logfddd=")
print.function(gev_p12_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12_ldda=function(x,t,v1,v2,v3,v4,v5){
	nx=length(x)
	ldd=matrix(0,5,5)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p12_logfdd)
	t1=t[,1]
	t2=t[,2]
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	ldd=deriv_copyldd(temp1,nx,dim=5)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12_lddda=function(x,t,v1,v2,v3,v4,v5){
	nx=length(x)
	lddd=array(0,c(5,5,5))

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p12_logfddd)
	t1=t[,1]
	t2=t[,2]
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	lddd=deriv_copylddd(temp1,nx,dim=5)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
