#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
#
f=function(x,v1,v2){(1/v2)*exp(-( ((x-v1)/v2) + exp(-((x-v1)/v2))))}
compare("d",extraDistr::dgumbel(1,2,3),f(1,2,3))
gumbel_fd=Deriv(f,c("v1","v2"),nderiv=1)
gumbel_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,v1,v2){exp(-(exp(-((x-v1)/v2))))}
compare("p",extraDistr::pgumbel(1,2,3),p(1,2,3))
gumbel_pd=Deriv(p,c("v1","v2"),nderiv=1)
gumbel_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2){-log(v2)-(x-v1)/v2-exp(-(x-v1)/v2)}
compare("l",extraDistr::dgumbel(1,2,3,log=TRUE),logf(1,2,3))
gumbel_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
gumbel_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("50c_gumbel_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gumbel_fd=")
print.function(gumbel_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gumbel_fdd=")
print.function(gumbel_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gumbel_pd=")
print.function(gumbel_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gumbel_pdd=")
print.function(gumbel_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gumbel_logfdd=")
print.function(gumbel_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gumbel_logfddd=")
print.function(gumbel_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gumbel_f1fa=function(x,v1,v2){
	vf=Vectorize(gumbel_fd,\"x\")
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
"gumbel_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_fdd,\"x\")
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
"gumbel_p1fa=function(x,v1,v2){
	vf=Vectorize(gumbel_pd,\"x\")
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
"gumbel_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_pdd,\"x\")
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
"gumbel_mu1fa=function(alpha,v1,v2){
	x=qgumbel((1-alpha),mu=v1,sigma=v2)
	vf=Vectorize(gumbel_pd,\"x\")
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
"gumbel_mu2fa=function(alpha,v1,v2){
	x=qgumbel((1-alpha),mu=v1,sigma=v2)
	nx=length(x)
	vf=Vectorize(gumbel_pdd,\"x\")
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
"gumbel_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_logfdd,\"x\")
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
"gumbel_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(gumbel_logfddd,\"x\")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
