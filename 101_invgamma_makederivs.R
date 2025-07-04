#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(actuar)
#
# v1 is shape (alpha)
# v2 is scale (theta)
f=function(x,v1,v2){((v2/x)^(v1))*exp(-(v2/x))*(1/x)*(1/gamma(v1))}
compare("d",actuar::dinvgamma(1,shape=2,scale=3),f(1,2,3))
invgamma_fd=Deriv(f,c("v1","v2"),nderiv=1)
invgamma_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
cat("  no easy cdf...involves regularized or incomplete gamma functions\n")
#
logf=function(x,v1,v2){v1*log(v2)-v1*log(x)-v2/x-log(x)-log(gamma(v1))}
compare("l",actuar::dinvgamma(1,shape=2,scale=3,log=TRUE),logf(1,2,3))
invgamma_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
invgamma_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("101c_invgamma_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("invgamma_fd=")
print.function(invgamma_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("invgamma_fdd=")
print.function(invgamma_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("invgamma_logfdd=")
print.function(invgamma_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("invgamma_logfddd=")
print.function(invgamma_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"invgamma_f1fa=function(x,v1,v2){
	vf=Vectorize(invgamma_fd,\"x\")
	f1=vf(x,v1,v2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"invgamma_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgamma_fdd,\"x\")
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"invgamma_p1fa=function(x,v1,v2){
###	vf=Vectorize(invgamma_pd,\"x\")
###	p1=vf(x,v1,v2)
###	return(p1)
###}\n"
###)
cat("############################################################\n")
#
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"invgamma_p2fa=function(x,v1,v2){
###	nx=length(x)
###	vf=Vectorize(invgamma_pdd,\"x\")
###	temp1=vf(x,v1,v2)
###	p2=deriv_copyfdd(temp1,nx,dim=2)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"invgamma_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgamma_logfdd,\"x\")
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"invgamma_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgamma_logfddd,\"x\")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
