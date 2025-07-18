#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
# v1=shape
# v2=scale
f=function(x,v1,v2){(v1/v2)*((x/v2)^(v1-1))*exp(-((x/v2)^v1))}
compare("d",dweibull(1,2,3),f(1,2,3))
weibull_fd=Deriv(f,c("v1","v2"),nderiv=1)
weibull_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,v1,v2){1-exp(-((x/v2)^v1))}
compare("p",pweibull(1,2,3),p(1,2,3))
weibull_pd=Deriv(p,c("v1","v2"),nderiv=1)
weibull_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2){log(v1)-log(v2)+(v1-1)*log(x/v2)-(x/v2)^v1}
compare("l",dweibull(1,2,3,log=TRUE),logf(1,2,3))
weibull_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
weibull_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("052c_weibull_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("weibull_fd=")
print.function(weibull_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("weibull_fdd=")
print.function(weibull_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("weibull_pd=")
print.function(weibull_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("weibull_pdd=")
print.function(weibull_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("weibull_logfdd=")
print.function(weibull_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("weibull_logfddd=")
print.function(weibull_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_f1fa=function(x,v1,v2){
	vf=Vectorize(weibull_fd,\"x\")
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
"weibull_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_fdd,\"x\")
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p1fa=function(x,v1,v2){
	vf=Vectorize(weibull_pd,\"x\")
	p1=vf(x,v1,v2)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_pdd,\"x\")
	temp1=vf(x,v1,v2)
	p2=deriv_copyfdd(temp1,nx,dim=2)
	return(p2)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_mu1fa=function(alpha,v1,v2){
	x=qweibull((1-alpha),shape=v1,scale=v2)
	vf=Vectorize(weibull_pd,\"x\")
	mu1=-vf(x,v1,v2)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_mu2fa=function(alpha,v1,v2){
	x=qweibull((1-alpha),shape=v1,scale=v2)
	nx=length(x)
	vf=Vectorize(weibull_pdd,\"x\")
	temp1=vf(x,v1,v2)
	mu2=-deriv_copyfdd(temp1,nx,dim=2)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_logfdd,\"x\")
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
"weibull_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(weibull_logfddd,\"x\")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
