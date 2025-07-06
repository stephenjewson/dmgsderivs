#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
#
f=function(x,t1,t2,v1,v2,v3,v4,v5){
	(exp(-v3-v4*t2))*((1+v5*((x-v1-v2*t1)*exp(-v3-v4*t2)))^(-1/v5-1))*exp(-((1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/v5)))
}
compare("d",extraDistr::dgev(1,4+5*2,exp(.6+.7*3),0.1),f(1,2,3,4,5,.6,.7,0.1))
gev_p12k3_fd=Deriv(f,c("v1","v2","v3","v4"),nderiv=1)
gev_p12k3_fdd=Deriv(f,c("v1","v2","v3","v4"),nderiv=2)
#
p=function(x,t1,t2,v1,v2,v3,v4,v5){
	exp(-((1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/v5)))
}
compare("p",extraDistr::pgev(1,4+5*2,exp(.6+.7*3),0.1),p(1,2,3,4,5,.6,.7,0.1))
gev_p12k3_pd=Deriv(p,c("v1","v2","v3","v4"),nderiv=1)
gev_p12k3_pdd=Deriv(p,c("v1","v2","v3","v4"),nderiv=2)
#
logf=function(x,t1,t2,v1,v2,v3,v4,v5){
	-v3-v4*t2-(1+1/v5)*log(1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))-(1+v5*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/v5)
}
compare("l",extraDistr::dgev(1,4+5*2,exp(.6+.7*3),0.1,log=TRUE),logf(1,2,3,4,5,.6,.7,0.1))
gev_p12k3_logfdd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=2)
gev_p12k3_logfddd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=3)
#
sink("074c_gev_p12k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p12k3_fd=")
print.function(gev_p12k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p12k3_fdd=")
print.function(gev_p12k3_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p12k3_pd=")
print.function(gev_p12k3_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p12k3_pdd=")
print.function(gev_p12k3_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p12k3_logfdd=")
print.function(gev_p12k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gev_p12k3_logfddd=")
print.function(gev_p12k3_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density for DMGS\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_f1fa=function(x,t01,t02,v1,v2,v3,v4,kshape){
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p12k3_fd,\"x\")
	f1=vf(x,t01,t02,v1,v2,v3,v4,kshape)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the density for WAIC\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_f1fw=function(x,t1,t2,v1,v2,v3,v4,kshape){
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p12k3_fd,c(\"x\",\"t1\",\"t2\"))
	f1=vf(x,t1,t2,v1,v2,v3,v4,kshape)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density for DMGS\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_f2fa=function(x,t01,t02,v1,v2,v3,v4,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_fdd,\"x\")
	temp1=vf(x,t01,t02,v1,v2,v3,v4,kshape)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density for WAIC\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_f2fw=function(x,t1,t2,v1,v2,v3,v4,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_fdd,c(\"x\",\"t1\",\"t2\"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4,kshape)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p12k3_p1fa=function(x,t01,t02,v1,v2,v3,v4,kshape){
###	kshape=movexiawayfromzero(kshape)
###	vf=Vectorize(gev_p12k3_pd,\"x\")
###	p1=vf(x,t01,t02,v1,v2,v3,v4,kshape)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p12k3_p2fa=function(x,t01,t02,v1,v2,v3,v4,kshape){
###	nx=length(x)
###
###	kshape=movexiawayfromzero(kshape)
###
###	vf=Vectorize(gev_p12k3_pdd,\"x\")
###	temp1=vf(x,t01,t02,v1,v2,v3,v4,kshape)
###	p2=deriv_copyfdd(temp1,nx,dim=4)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_mu1fa=function(alpha,t01,t02,v1,v2,v3,v4,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=kshape)
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p12k3_pd,\"x\")
	mu1=-vf(x,t01,t02,v1,v2,v3,v4,kshape)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_mu2fa=function(alpha,t01,t02,v1,v2,v3,v4,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=kshape)
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_pdd,\"x\")
	temp1=vf(x,t01,t02,v1,v2,v3,v4,kshape)
	mu2=-deriv_copyfdd(temp1,nx,dim=4)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_ldda=function(x,t1,t2,v1,v2,v3,v4,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p12k3_logfdd,c(\"x\",\"t1\",\"t2\"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4,kshape)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p12k3_lddda=function(x,t1,t2,v1,v2,v3,v4,kshape){
	nx=length(x)
	vf=Vectorize(gev_p12k3_logfddd,c(\"x\",\"t1\",\"t2\"))

	kshape=movexiawayfromzero(kshape)

	temp1=vf(x,t1,t2,v1,v2,v3,v4,kshape)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
