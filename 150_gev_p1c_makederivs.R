#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
f=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	(1/v5)*((1+v6*((x-v1-v2*ta-v3*tb-v4*tc)/v5))^(-1/v6-1))*
		exp(-((1+v6*(x-v1-v2*ta-v3*tb-v4*tc)/v5)^(-1/v6)))
}
compare("d",extraDistr::dgev(1,3+4*2,5,0.1),f(1,2,0,0,3,4,0,0,5,0.1))
gev_p1c_fd=Deriv(f,c("v1","v2","v3","v4","v5","v6"),nderiv=1)
gev_p1c_fdd=Deriv(f,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
#
p=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	exp(-((1+v6*(x-v1-v2*ta-v3*tb-v4*tc)/v5)^(-1/v6)))
	}
compare("p",extraDistr::pgev(1,3+4*2,5,0.1),p(1,2,0,0,3,4,0,0,5,0.1))
gev_p1c_pd=Deriv(p,c("v1","v2","v3","v4","v5","v6"),nderiv=1)
gev_p1c_pdd=Deriv(p,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
#
logf=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	-log(v5)-(1+1/v6)*log(1+v6*(x-v1-v2*ta-v3*tb-v4*tc)/v5)-(1+v6*(x-v1-v2*ta-v3*tb-v4*tc)/v5)^(-1/v6)
}
compare("l",extraDistr::dgev(1,3+4*2,5,0.1,log=TRUE),logf(1,2,0,0,3,4,0,0,5,0.1))
gev_p1c_logfdd=Deriv(logf,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
gev_p1c_logfddd=Deriv(logf,c("v1","v2","v3","v4","v5","v6"),nderiv=3)
#
sink("150c_gev_p1c_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p1c_fd=")
print.function(gev_p1c_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1c_fdd=")
print.function(gev_p1c_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p1c_pd=")
print.function(gev_p1c_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1c_pdd=")
print.function(gev_p1c_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1c_logfdd=")
print.function(gev_p1c_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gev_p1c_logfddd=")
print.function(gev_p1c_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density for DMGS\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_f1fa=function(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fd,\"x\")
	f1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the density for WAIC\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_f1fw=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fd,c(\"x\",\"ta\",\"tb\",\"tc\"))
	f1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density for DMGS\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_f2fa=function(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fdd,\"x\")
	temp1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	f2=deriv_copyfdd(temp1,nx,dim=6)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density for WAIC\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_f2fw=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fdd,c(\"x\",\"ta\",\"tb\",\"tc\"))
	temp1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	f2=deriv_copyfdd(temp1,nx,dim=6)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1c_p1fa=function(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
###
###	v6=movexiawayfromzero(v6)
###
###	vf=Vectorize(gev_p1c_pd,\"x\")
###	p1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1c_p2fa=function(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
###	nx=length(x)
###
###	v6=movexiawayfromzero(v6)
###
###	vf=Vectorize(gev_p1c_pdd,\"x\")
###	temp1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
###	p2=deriv_copyfdd(temp1,nx,dim=6)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_mu1fa=function(alpha,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b+v4*t0c,sigma=v5,xi=v6)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_pd,\"x\")
	mu1=-vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_mu2fa=function(alpha,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b+v4*t0c,sigma=v5,xi=v6)
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_pdd,\"x\")
	temp1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	mu2=-deriv_copyfdd(temp1,nx,dim=6)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_ldda=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_logfdd,c(\"x\",\"ta\",\"tb\",\"tc\"))
	temp1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	ldd=deriv_copyldd(temp1,nx,dim=6)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1c_lddda=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_logfddd,c(\"x\",\"ta\",\"tb\",\"tc\"))
	temp1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	lddd=deriv_copylddd(temp1,nx,dim=6)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
