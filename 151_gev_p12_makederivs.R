#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
f=function(x,t1,t2,v1,v2,v3,v4,v5){
	(exp(-v3-v4*t2))*((1+v5*((x-v1-v2*t1)*exp(-v3-v4*t2)))^(-1/v5-1))*exp(-((1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/v5)))
}
compare("d",extraDistr::dgev(1,4+5*2,exp(.6+.7*3),0.1),f(1,2,3,4,5,.6,.7,0.1))
gev_p12_fd=Deriv(f,c("v1","v2","v3","v4","v5"),nderiv=1)
gev_p12_fdd=Deriv(f,c("v1","v2","v3","v4","v5"),nderiv=2)
#
p=function(x,t1,t2,v1,v2,v3,v4,v5){
	exp(-((1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/v5)))
}
compare("p",extraDistr::pgev(1,4+5*2,exp(.6+.7*3),0.1),p(1,2,3,4,5,.6,.7,0.1))
gev_p12_pd=Deriv(p,c("v1","v2","v3","v4","v5"),nderiv=1)
gev_p12_pdd=Deriv(p,c("v1","v2","v3","v4","v5"),nderiv=2)
#
logf=function(x,t1,t2,v1,v2,v3,v4,v5){
	-v3-v4*t2-(1+1/v5)*log(1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))-(1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/v5)
}
compare("l",extraDistr::dgev(1,4+5*2,exp(.6+.7*3),0.1,log=TRUE),logf(1,2,3,4,5,.6,.7,0.1))
gev_p12_logfdd=Deriv(logf,c("v1","v2","v3","v4","v5"),nderiv=2)
gev_p12_logfddd=Deriv(logf,c("v1","v2","v3","v4","v5"),nderiv=3)
#
sink("151c_gev_p12_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_fd=")
print.function(gev_p12_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_fdd=")
print.function(gev_p12_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_pd=")
print.function(gev_p12_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_pdd=")
print.function(gev_p12_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_logfdd=")
print.function(gev_p12_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gev_p12_logfddd=")
print.function(gev_p12_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12_f1fa=function(x,t01,t02,v1,v2,v3,v4,v5){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_fd,\"x\")
	f1=vf(x,t01,t02,v1,v2,v3,v4,v5)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12_f2fa=function(x,t01,t02,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_fdd,\"x\")
	temp1=vf(x,t01,t02,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p12_p1fa=function(x,t01,t02,v1,v2,v3,v4,v5){
###
###	v3=movexiawayfromzero(v3)
###
###	vf=Vectorize(gev_p12_pd,\"x\")
###	p1=vf(x,t01,t02,v1,v2,v3,v4,v5)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p12_p2fa=function(x,t01,t02,v1,v2,v3,v4,v5){
###	nx=length(x)
###
###	v3=movexiawayfromzero(v3)
###
###	vf=Vectorize(gev_p12_pdd,\"x\")
###	temp1=vf(x,t01,t02,v1,v2,v3,v4,v5)
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
"gev_p12_mu1fa=function(alpha,t01,t02,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=v5)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_pd,\"x\")
	mu1=-vf(x,t01,t02,v1,v2,v3,v4,v5)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12_mu2fa=function(alpha,t01,t02,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=v5)
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_pdd,\"x\")
	temp1=vf(x,t01,t02,v1,v2,v3,v4,v5)
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
"gev_p12_ldda=function(x,t1,t2,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_logfdd,\"x\")
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
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
"gev_p12_lddda=function(x,t1,t2,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_logfddd,\"x\")
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	lddd=deriv_copylddd(temp1,nx,dim=5)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
