#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1,v2){(1/sqrt(2*pi))*(1/v2)*exp(-(x-v1)^2/(2*v2*v2))}
compare("d",dnorm(1,2,3),f(1,2,3))
norm_fd=Deriv(f,c("v1","v2"),nderiv=1)
norm_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
cat("  no cdf\n")
#
logf=function(x,v1,v2){-0.5*log(2*pi)-log(v2)-(x-v1)^2/(2*v2*v2)}
compare("l",dnorm(1,2,3,log=TRUE),logf(1,2,3))
norm_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
norm_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("30c_norm_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_fd=")
print.function(norm_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_fdd=")
print.function(norm_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_logfdd=")
print.function(norm_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_logfddd=")
print.function(norm_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"norm_f1fa=function(x,v1,v2){
	vf=Vectorize(norm_fd)
	f1=vf(x,v1,v2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"norm_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_logfdd)
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
"norm_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(norm_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
