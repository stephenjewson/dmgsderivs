#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2){exp(-v1-v2*t)*exp(-x*exp(-v1-v2*t))}
compare("d",dexp(1,exp(-3-4*2)),f(1,2,3,4))
exp_p1_fd=Deriv(f,c("v1","v2"),nderiv=1)
exp_p1_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,t,v1,v2){1-exp(-x*exp(-v1-v2*t))}
compare("p",pexp(1,exp(-3-4*2)),p(1,2,3,4))
exp_p1_pd=Deriv(p,c("v1","v2"),nderiv=1)
exp_p1_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
#logf=function(x,t,v1,v2){log(v1)-v1*x}
logf=function(x,t,v1,v2){-v1-v2*t-x*exp(-v1-v2*t)}
compare("l",dexp(1,exp(-3-4*2),log=TRUE),logf(1,2,3,4))
exp_p1_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
exp_p1_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("055c_exp_p1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_fd=")
print.function(exp_p1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_fdd=")
print.function(exp_p1_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_pd=")
print.function(exp_p1_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_pdd=")
print.function(exp_p1_pdd)
cat("######################################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_logfdd=")
print.function(exp_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("exp_p1_logfddd=")
print.function(exp_p1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_f1fa=function(x,t,v1,v2){
	vf=Vectorize(exp_p1_fd,\"x\")
	f1=vf(x,t,v1,v2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_f2fa=function(x,t,v1,v2){
	nx=length(x)
	vf=Vectorize(exp_p1_fdd,\"x\")
	temp1=vf(x,t,v1,v2)
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
"exp_p1_p1fa=function(x,t,v1,v2){
	vf=Vectorize(exp_p1_pd,\"x\")
	p1=vf(x,t,v1,v2)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_p2fa=function(x,t,v1,v2){
	nx=length(x)
	p2=array(0,c(2,2,nx))
	vf=Vectorize(exp_p1_pdd,\"x\")
	temp1=vf(x,t,v1,v2)
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
"exp_p1_mu1fa=function(alpha,t,v1,v2){
	x=qexp((1-alpha),rate=exp(-v1-v2*t))
	vf=Vectorize(exp_p1_pd,\"x\")
	mu1=-vf(x,t,v1,v2)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_mu2fa=function(alpha,t,v1,v2){
	x=qexp((1-alpha),rate=exp(-v1-v2*t))
	nalpha=length(alpha)
	vf=Vectorize(exp_p1_pdd,\"x\")
	temp1=vf(x,t,v1,v2)
	mu2=-deriv_copyfdd(temp1,nalpha,dim=2)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1_ldda=function(x,t,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(exp_p1_logfdd,\"x\")
	temp1=vf(x,t,v1,v2)
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
"exp_p1_lddda=function(x,t,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(exp_p1_logfddd,\"x\")
	temp1=vf(x,t,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))

