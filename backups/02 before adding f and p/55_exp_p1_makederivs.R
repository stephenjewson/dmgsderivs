#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2){-v1-v2*t-x*exp(-v1-v2*t)}
exp_p1_logfdd=Deriv(f,c("v1","v2"),nderiv=2)
exp_p1_logfddd=Deriv(f,c("v1","v2"),nderiv=3)
#
sink("55c_exp_p1_derivs.R")
#
cat("######################################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_logfdd=")
print.function(exp_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_logfddd=")
print.function(exp_p1_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_ldda=function(x,t,v1,v2){
	nx=length(x)
	vf=Vectorize(exp_p1_logfdd)
	temp1=vf(x,t,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_lddda=function(x,t,v1,v2){
	nx=length(x)
	vf=Vectorize(exp_p1_logfddd)
	temp1=vf(x,t,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

