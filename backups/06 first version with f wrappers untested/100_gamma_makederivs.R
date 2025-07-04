#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# v1 is shape (alpha)
# v2 is scale (sigma)
f=function(x,v1,v2){(v2^(-v1))*(1/gamma(v1))*(x^(v1-1))*exp(-(x/v2))}
compare("d",dgamma(1,shape=2,scale=3),f(1,2,3))
gamma_fd=Deriv(f,c("v1","v2"),nderiv=1)
gamma_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
cat("  no cdf...involves lower incomplete gamma function\n")
#
logf=function(x,v1,v2){-v1*log(v2)-log(gamma(v1))+(v1-1)*log(x)-x/v2}
compare("l",dgamma(1,shape=2,scale=3,log=TRUE),logf(1,2,3))
gamma_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
gamma_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("100c_gamma_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gamma_fd=")
print.function(gamma_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gamma_fdd=")
print.function(gamma_fdd)
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
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gamma_f1fa=function(x,v1,v2){
	nx=length(x)
	f1=matrix(0,2,nx)
	vf=Vectorize(gamma_fd)
	temp1=vf(x,v1,v2)
	f1=deriv_copyldd(temp1,nx,dim=2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gamma_f2fa=function(x,v1,v2){
	nx=length(x)
	f2=array(0,c(2,2,nx))
	vf=Vectorize(gamma_fdd)
	temp1=vf(x,v1,v2)
	f2f=deriv_copyldd(temp1,nx,dim=2)
	return(f2f)
}\n"
)
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
