#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
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
cat("#' @inheritParams manf\n")
cat("invgamma_fd=")
print.function(invgamma_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("ingvamma_fdd=")
print.function(invgamma_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("invgamma_logfdd=")
print.function(invgamma_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("invgamma_logfddd=")
print.function(invgamma_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"invgamma_ldda=function(x,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(invgamma_logfdd)
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
"invgamma_lddda=function(x,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(invgamma_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
