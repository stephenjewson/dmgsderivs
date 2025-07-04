#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(gnorm)
#
f=function(x,v1,v2,v3){(1/(2*v2))*(1/gamma(1/v3))*v3*exp(-(abs(x-v1)/v2)^v3)}
compare("d",gnorm::dgnorm(1,2,3,4),f(1,2,3,4))
gnorm_k3_fd=Deriv(f,c("v1","v2"),nderiv=1)
gnorm_k3_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
cat("  no cdf\n")
#
logf=function(x,v1,v2,v3){log(v3)-((abs(x-v1))/v2)^v3-log(2)-log(v2)-log(gamma(1/v3))}
compare("l",gnorm::dgnorm(1,2,3,4,log=TRUE),logf(1,2,3,4))
gnorm_k3_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
gnorm_k3_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("32c_gnorm_k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gnorm_k3_fd=")
print.function(gnorm_k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gnorm_k3_fdd=")
print.function(gnorm_k3_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gnorm_k3_logfdd=")
print.function(gnorm_k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gnorm_k3_logfddd=")
print.function(gnorm_k3_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gnorm_k3_f1fa=function(x,v1,v2,kbeta){
	vf=Vectorize(gnorm_k3_fd)
	f1=vf(x,v1,v2,kbeta)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gnorm_k3_f2fa=function(x,v1,v2,kbeta){
	nx=length(x)
	vf=Vectorize(gnorm_k3_fdd)
	temp1=vf(x,v1,v2,kbeta)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gnorm_k3_ldda=function(x,v1,v2,kbeta){
	nx=length(x)
	vf=Vectorize(gnorm_k3_logfdd)
	temp1=vf(x,v1,v2,kbeta)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gnorm_k3_lddda=function(x,v1,v2,kbeta){
	nx=length(x)
	vf=Vectorize(gnorm_k3_logfddd)
	temp1=vf(x,v1,v2,kbeta)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

