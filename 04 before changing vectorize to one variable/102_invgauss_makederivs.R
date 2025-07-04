#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(actuar)
#
# note that I'm using shape (=1/dispersion)
f=function(x,v1,v2){((v2/(2*pi*x*x*x))^(0.5))*exp(-(v2*(x-v1)^2/(2*v1*v1*x)))}
compare("d",actuar::dinvgauss(1,2,shape=3),f(1,2,3))
invgauss_fd=Deriv(f,c("v1","v2"),nderiv=1)
invgauss_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
cat("  no cdf (well, involves erf)\n")
#
logf=function(x,v1,v2){0.5*log(v2)-0.5*log(2*pi)-1.5*log(x)-v2*(x-v1)^2/(2*v1*v1*x)}
compare("l",actuar::dinvgauss(1,2,shape=3,log=TRUE),logf(1,2,3))
invgauss_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
invgauss_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("102c_invgauss_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("invgauss_fd=")
print.function(invgauss_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("invgauss_fdd=")
print.function(invgauss_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("invgauss_logfdd=")
print.function(invgauss_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("invgauss_logfddd=")
print.function(invgauss_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"invgauss_f1fa=function(x,v1,v2){
	vf=Vectorize(invgauss_fd,"x")
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
"invgauss_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgauss_fdd,"x")
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
###"invgauss_p1fa=function(x,v1,v2){
###	vf=Vectorize(invgauss_pd,"x")
###	p1=vf(x,v1,v2)
###	return(p1)
###}\n"
###)
cat("############################################################\n")
#
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"invgauss_p2fa=function(x,v1,v2){
###	nx=length(x)
###	vf=Vectorize(invgauss_pdd,"x")
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
"invgauss_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgauss_logfdd,"x")
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
"invgauss_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(invgauss_logfddd,"x")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/makederivatives/',sep=""))
