#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2,v3){(1/(pi*v3))*1/(1+(x-v1-v2*t)^2/(v3*v3))}
compare("d",dcauchy(1,3+4*2,5),f(1,2,3,4,5))
cauchy_p1_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
cauchy_p1_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
cat("  no cdf (well, contains arctan)\n")
#
logf=function(x,t,v1,v2,v3){-log(pi)-log(v3)-log(1+((x-v1-v2*t)/v3)^2)}
compare("l",dcauchy(1,3+4*2,5,log=TRUE),logf(1,2,3,4,5))
cauchy_p1_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
cauchy_p1_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("64c_cauchy_p1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("cauchy_p1_fd=")
print.function(cauchy_p1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("cauchy_p1_fdd=")
print.function(cauchy_p1_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("cauchy_p1_logfdd=")
print.function(cauchy_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("cauchy_p1_logfddd=")
print.function(cauchy_p1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"cauchy_p1_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(cauchy_p1_fd,\"x\")
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"cauchy_p1_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_fdd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"cauchy_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_logfdd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"cauchy_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(cauchy_p1_logfddd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
