#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1,v2){-v1*log(v2)-log(gamma(v1))+(v1-1)*log(x)-x/v2}
gamma_logfdd=Deriv(f,c("v1","v2"),nderiv=2)
gamma_logfddd=Deriv(f,c("v1","v2"),nderiv=3)
#
sink("100c_gamma_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gamma_logfdd=")
print.function(gamma_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gamma_logfddd=")
print.function(gamma_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gamma_ldda=function(x,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(gamma_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gamma_lddda=function(x,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(gamma_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
