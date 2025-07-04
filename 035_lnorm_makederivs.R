#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
#
f=function(x,v1,v2){(1/sqrt(2*pi))*(1/v2)*(1/x)*exp(-(log(x)-v1)^2/(2*v2*v2))}
compare("d",dlnorm(1,2,3),f(1,2,3))
lnorm_fd=Deriv(f,c("v1","v2"),nderiv=1)
lnorm_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,v1,v2){pnorm((log(x)-v1)/v2)}
compare("p",plnorm(1,2,3),p(1,2,3))
lnorm_pd=Deriv(p,c("v1","v2"),nderiv=1)
lnorm_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2){-0.5*log(2*pi)-log(v2)-log(x)-(log(x)-v1)^2/(2*v2*v2)}
compare("l",dlnorm(1,2,3,log=TRUE),logf(1,2,3))
lnorm_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
lnorm_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("035c_lnorm_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("lnorm_fd=")
print.function(lnorm_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("lnorm_fdd=")
print.function(lnorm_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("lnorm_pd=")
print.function(lnorm_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("lnorm_pdd=")
print.function(lnorm_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("lnorm_logfdd=")
print.function(lnorm_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("lnorm_logfddd=")
print.function(lnorm_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"lnorm_f1fa=function(x,v1,v2){
	vf=Vectorize(lnorm_fd,\"x\")
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
"lnorm_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_fdd,\"x\")
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
"lnorm_p1fa=function(x,v1,v2){
	vf=Vectorize(lnorm_pd,\"x\")
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
"lnorm_p2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_pdd,\"x\")
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
"lnorm_mu1fa=function(alpha,v1,v2){
	x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
	vf=Vectorize(lnorm_pd,\"x\")
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
"lnorm_mu2fa=function(alpha,v1,v2){
	x=qlnorm((1-alpha),meanlog=v1,sdlog=v2)
	nx=length(x)
	vf=Vectorize(lnorm_pdd,\"x\")
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
"lnorm_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_logfdd,\"x\")
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
"lnorm_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(lnorm_logfddd,\"x\")
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))

