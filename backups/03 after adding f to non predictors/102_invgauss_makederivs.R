#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(actuar)
#
# note that I'm using shape (=1/dispersion)
f=function(x,v1,v2){((v2/(2*pi*x*x*x))^(0.5))*exp(-(v2*(x-v1)^2/(2*v1*v1*x)))}
cat(actuar::dinvgauss(1,2,shape=3),":",f(1,2,3),"\n")
invgauss_fd=Deriv(f,c("v1","v2"),nderiv=1)
invgauss_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2){0.5*log(v2)-0.5*log(2*pi)-1.5*log(x)-v2*(x-v1)^2/(2*v1*v1*x)}
invgauss_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
invgauss_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("102c_invgauss_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("invgauss_fd=")
print.function(invgauss_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("invgauss_fdd=")
print.function(invgauss_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("invgauss_logfdd=")
print.function(invgauss_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("invgauss_logfddd=")
print.function(invgauss_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"invgauss_ldda=function(x,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(invgauss_logfdd)
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
"invgauss_lddda=function(x,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(invgauss_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
