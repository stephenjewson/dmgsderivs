#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
#
f=function(x,t,v1,v2,v3,v4){(1/v3)*((1+v4*((x-v1-v2*t)/v3))^(-1/v4-1))*exp(-((1+v4*(x-v1-v2*t)/v3)^(-1/v4)))}
compare("d",extraDistr::dgev(1,3+4*2,5,0.1),f(1,2,3,4,5,0.1))
gev_p1k3_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
gev_p1k3_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
p=function(x,t,v1,v2,v3,v4){exp(-((1+v4*(x-v1-v2*t)/v3)^(-1/v4)))}
compare("p",extraDistr::pgev(1,3+4*2,5,0.1),p(1,2,3,4,5,0.1))
gev_p1k3_pd=Deriv(p,c("v1","v2","v3"),nderiv=1)
gev_p1k3_pdd=Deriv(p,c("v1","v2","v3"),nderiv=2)
#
logf=function(x,t,v1,v2,v3,v4){-log(v3)-(1+1/v4)*log(1+v4*(x-v1-v2*t)/v3)-(1+v4*(x-v1-v2*t)/v3)^(-1/v4)}
compare("l",extraDistr::dgev(1,3+4*2,5,0.1,log=TRUE),logf(1,2,3,4,5,0.1))
gev_p1k3_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
gev_p1k3_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("74c_gev_p1k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p1k3_fd=")
print.function(gev_p1k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1k3_fdd=")
print.function(gev_p1k3_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gev_p1k3_pd=")
print.function(gev_p1k3_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1k3_pdd=")
print.function(gev_p1k3_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gev_p1k3_logfdd=")
print.function(gev_p1k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gev_p1k3_logfddd=")
print.function(gev_p1k3_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1k3_f1fa=function(x,t,v1,v2,v3,kshape){
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p1k3_fd,"x")
	f1=vf(x,t,v1,v2,v3,kshape)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1k3_f2fa=function(x,t,v1,v2,v3,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p1k3_fdd,"x")
	temp1=vf(x,t,v1,v2,v3,kshape)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1k3_p1fa=function(x,t,v1,v2,v3,kshape){
###	kshape=movexiawayfromzero(kshape)
###	vf=Vectorize(gev_p1k3_pd,"x")
###	p1=vf(x,t,v1,v2,v3,kshape)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1k3_p2fa=function(x,t,v1,v2,v3,kshape){
###	nx=length(x)
###
###	kshape=movexiawayfromzero(kshape)
###
###	vf=Vectorize(gev_p1k3_pdd,"x")
###	temp1=vf(x,t,v1,v2,v3,kshape)
###	p2=deriv_copyfdd(temp1,nx,dim=3)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1k3_mu1fa=function(alpha,t,v1,v2,v3,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=kshape)
	kshape=movexiawayfromzero(kshape)
	vf=Vectorize(gev_p1k3_pd,"x")
	mu1=-vf(x,t,v1,v2,v3,kshape)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1k3_mu2fa=function(alpha,t,v1,v2,v3,kshape){
	x=extraDistr::qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=kshape)
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p1k3_pdd,"x")
	temp1=vf(x,t,v1,v2,v3,kshape)
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
"gev_p1k3_ldda=function(x,t,v1,v2,v3,kshape){
	nx=length(x)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_p1k3_logfdd,"x")
	temp1=vf(x,t,v1,v2,v3,kshape)
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
"gev_p1k3_lddda=function(x,t,v1,v2,v3,kshape){
	nx=length(x)
	vf=Vectorize(gev_p1k3_logfddd,"x")

	kshape=movexiawayfromzero(kshape)

	temp1=vf(x,t,v1,v2,v3,kshape)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/makederivatives/',sep=""))
