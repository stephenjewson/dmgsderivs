#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
#
f=function(x,t,v1,v2,v3){(1/v3)*exp(-( ((x-v1-v2*t)/v3) + exp(-((x-v1-v2*t)/v3))))}
compare("d",extraDistr::dgumbel(1,3+4*2,5),f(1,2,3,4,5))
gumbel_p1_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
gumbel_p1_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
p=function(x,t,v1,v2,v3){exp(-(exp(-((x-v1-v2*t)/v3))))}
compare("p",extraDistr::pgumbel(1,3+4*2,5),p(1,2,3,4,5))
gumbel_p1_pd=Deriv(p,c("v1","v2","v3"),nderiv=1)
gumbel_p1_pdd=Deriv(p,c("v1","v2","v3"),nderiv=2)
#
logf=function(x,t,v1,v2,v3){-log(v3)-(x-v1-v2*t)/v3-exp(-(x-v1-v2*t)/v3)}
compare("l",extraDistr::dgumbel(1,3+4*2,5,log=TRUE),logf(1,2,3,4,5))
gumbel_p1_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
gumbel_p1_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("70c_gumbel_p1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gumbel_p1_fd=")
print.function(gumbel_p1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gumbel_p1_fdd=")
print.function(gumbel_p1_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("gumbel_p1_pd=")
print.function(gumbel_p1_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gumbel_p1_pdd=")
print.function(gumbel_p1_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("gumbel_p1_logfdd=")
print.function(gumbel_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("gumbel_p1_logfddd=")
print.function(gumbel_p1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"gumbel_p1_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(gumbel_p1_fd,\"x\")
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
"gumbel_p1_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_fdd,\"x\")
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
"gumbel_p1_p1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(gumbel_p1_pd,\"x\")
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
"gumbel_p1_p2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_pdd,\"x\")
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
"gumbel_p1_mu1fa=function(alpha,t,v1,v2,v3){
	x=qgumbel((1-alpha),mu=v1+v2*t,sigma=v3)
	vf=Vectorize(gumbel_p1_pd,\"x\")
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
"gumbel_p1_mu2fa=function(alpha,t,v1,v2,v3){
	x=qgumbel((1-alpha),mu=v1+v2*t,sigma=v3)
	nx=length(x)
	vf=Vectorize(gumbel_p1_pdd,\"x\")
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
"gumbel_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_logfdd,\"x\")
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
"gumbel_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(gumbel_p1_logfddd,\"x\")
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
