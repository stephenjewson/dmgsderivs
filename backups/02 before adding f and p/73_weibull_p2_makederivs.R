#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2,v3){log(v1)-v2-t*v3+(v1-1)*(log(x)-v2-t*v3)-(x*exp(-v2-t*v3))^v1}
weibull_p2_logfdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
weibull_p2_logfddd=Deriv(f,c("v1","v2","v3"),nderiv=3)
#
sink("73c_weibull_p2_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_logfdd=")
print.function(weibull_p2_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_logfddd=")
print.function(weibull_p2_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	ldd=matrix(0,3,3)
	vf=Vectorize(weibull_p2_logfdd)
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	lddd=array(0,c(3,3,3))
	vf=Vectorize(weibull_p2_logfddd)
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
