#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
f=function(x,ta,tb,v1,v2,v3,v4,v5){
	(1/v4)*((1+v5*((x-v1-v2*ta-v3*tb)/v4))^(-1/v5-1))*exp(-((1+v5*(x-v1-v2*ta-v3*tb)/v4)^(-1/v5)))
}
compare("d",extraDistr::dgev(1,3+4*2,5,0.1),f(1,2,0,3,4,0,5,0.1))
gev_p1b_fd=Deriv(f,c("v1","v2","v3","v4","v5"),nderiv=1)
gev_p1b_fdd=Deriv(f,c("v1","v2","v3","v4","v5"),nderiv=2)
#
p=function(x,ta,tb,v1,v2,v3,v4,v5){
	exp(-((1+v5*(x-v1-v2*ta-v3*tb)/v4)^(-1/v5)))
	}
compare("p",extraDistr::pgev(1,3+4*2,5,0.1),p(1,2,0,3,4,0,5,0.1))
gev_p1b_pd=Deriv(p,c("v1","v2","v3","v4","v5"),nderiv=1)
gev_p1b_pdd=Deriv(p,c("v1","v2","v3","v4","v5"),nderiv=2)
#
logf=function(x,ta,tb,v1,v2,v3,v4,v5){
	-log(v4)-(1+1/v5)*log(1+v5*(x-v1-v2*ta-v3*tb)/v4)-(1+v5*(x-v1-v2*ta-v3*tb)/v4)^(-1/v5)
}
compare("l",extraDistr::dgev(1,3+4*2,5,0.1,log=TRUE),logf(1,2,0,3,4,0,5,0.1))
gev_p1b_logfdd=Deriv(logf,c("v1","v2","v3","v4","v5"),nderiv=2)
gev_p1b_logfddd=Deriv(logf,c("v1","v2","v3","v4","v5"),nderiv=3)
#
sink("150c_gev_p1b_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p1b_fd=")
print.function(gev_p1b_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1b_fdd=")
print.function(gev_p1b_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p1b_pd=")
print.function(gev_p1b_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1b_pdd=")
print.function(gev_p1b_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1b_logfdd=")
print.function(gev_p1b_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gev_p1b_logfddd=")
print.function(gev_p1b_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density for DMGS\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_f1fa=function(x,t0a,t0b,v1,v2,v3,v4,v5){

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fd,\"x\")
	f1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the density for WAIC\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_f1fw=function(x,ta,tb,v1,v2,v3,v4,v5){

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fd,c(\"x\",\"ta\",\"tb\"))
	f1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density for DMGS\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_f2fa=function(x,t0a,t0b,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fdd,\"x\")
	temp1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density for WAIC\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_f2fw=function(x,ta,tb,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fdd,c(\"x\",\"ta\",\"tb\"))
	temp1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1b_p1fa=function(x,t0a,t0b,v1,v2,v3,v4,v5){
###
###	v5=movexiawayfromzero(v5)
###
###	vf=Vectorize(gev_p1b_pd,\"x\")
###	p1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1b_p2fa=function(x,t0a,t0b,v1,v2,v3,v4,v5){
###	nx=length(x)
###
###	v5=movexiawayfromzero(v5)
###
###	vf=Vectorize(gev_p1b_pdd,\"x\")
###	temp1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
###	p2=deriv_copyfdd(temp1,nx,dim=5)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_mu1fa=function(alpha,t0a,t0b,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b,sigma=v4,xi=v5)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_pd,\"x\")
	mu1=-vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_mu2fa=function(alpha,t0a,t0b,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b,sigma=v4,xi=v5)
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_pdd,\"x\")
	temp1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	mu2=-deriv_copyfdd(temp1,nx,dim=5)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_ldda=function(x,ta,tb,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_logfdd,c(\"x\",\"ta\",\"tb\"))
	temp1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	ldd=deriv_copyldd(temp1,nx,dim=5)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1b_lddda=function(x,ta,tb,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_logfddd,c(\"x\",\"ta\",\"tb\"))
	temp1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	lddd=deriv_copylddd(temp1,nx,dim=5)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
