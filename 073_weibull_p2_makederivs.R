#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
# v1=shape
# v2=scale
f=function(x,t,v1,v2,v3){(v1/exp(v2+v3*t))*((x/exp(v2+v3*t))^(v1-1))*exp(-((x/exp(v2+v3*t))^v1))}
compare("d",dweibull(1,3,exp(.4+.5*.2)),f(1,.2,3,.4,.5))
weibull_p2_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
weibull_p2_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
p=function(x,t,v1,v2,v3){1-exp(-((x*exp(-v2-v3*t))^v1))}
compare("p",pweibull(1,3,exp(.4+.5*.2)),p(1,.2,3,.4,.5))
weibull_p2_pd=Deriv(p,c("v1","v2","v3"),nderiv=1)
weibull_p2_pdd=Deriv(p,c("v1","v2","v3"),nderiv=2)
#
logf=function(x,t,v1,v2,v3){log(v1)-v2-v3*t+(v1-1)*log(x)-(v1-1)*(v2+v3*t)-(x/exp(v2+v3*t))^v1}
compare("l",dweibull(1,3,exp(.4+.5*.2),log=TRUE),logf(1,.2,3,.4,.5))
weibull_p2_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
weibull_p2_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("073c_weibull_p2_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_fd=")
print.function(weibull_p2_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_fdd=")
print.function(weibull_p2_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_pd=")
print.function(weibull_p2_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_pdd=")
print.function(weibull_p2_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_logfdd=")
print.function(weibull_p2_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("weibull_p2_logfddd=")
print.function(weibull_p2_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(weibull_p2_fd,\"x\")
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
"weibull_p2_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_fdd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_p1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(weibull_p2_pd,\"x\")
	p1=vf(x,t,v1,v2,v3)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_p2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_pdd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_mu1fa=function(alpha,t,v1,v2,v3){
	x=qweibull((1-alpha),shape=v1,scale=exp(v2+v3*t))
	vf=Vectorize(weibull_p2_pd,\"x\")
	mu1=-vf(x,t,v1,v2,v3)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_mu2fa=function(alpha,t,v1,v2,v3){
	x=qweibull((1-alpha),shape=v1,scale=exp(v2+v3*t))
	nx=length(x)
	vf=Vectorize(weibull_p2_pdd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_p2_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_logfdd,\"x\")
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
"weibull_p2_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(weibull_p2_logfddd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
