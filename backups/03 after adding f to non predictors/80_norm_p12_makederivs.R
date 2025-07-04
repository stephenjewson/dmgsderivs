#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# mu=v1 + b v2
# sd=exp(v3+d v4)
logf=function(x,t1,t2,v1,v2,v3,v4){-0.5*log(2*pi)-v3-t2*v4-(x-v1-t1*v2)^2/(2*exp(2*v3+2*t2*v4))}
norm_p12_logfdd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=2)
norm_p12_logfddd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=3)
#
sink("80c_norm_p12_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_logfdd=")
print.function(norm_p12_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_logfddd=")
print.function(norm_p12_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_ldda=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	ldd=matrix(0,4,4)
	vf=Vectorize(norm_p12_logfdd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_lddda=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	lddd=array(0,c(4,4,4))
	vf=Vectorize(norm_p12_logfddd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
