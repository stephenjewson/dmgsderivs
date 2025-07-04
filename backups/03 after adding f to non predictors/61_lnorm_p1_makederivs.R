#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
logf=function(x,t,v1,v2,v3){-0.5*log(2*pi)-log(v3)-log(x)-(log(x)-v1-t*v2)^2/(2*v3*v3)}
lnorm_p1_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
lnorm_p1_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("61c_lnorm_p1_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lnorm_p1_logfdd=")
print.function(lnorm_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lnorm_p1_logfddd=")
print.function(lnorm_p1_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lnorm_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	ldd=matrix(0,3,3)
	vf=Vectorize(lnorm_p1_logfdd)
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
"lnorm_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	lddd=array(0,c(3,3,3))
	vf=Vectorize(lnorm_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

