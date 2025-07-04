#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t1,t2,v1,v2,v3,v4){(1/sqrt(2*pi))*(1/exp(v3+v4*t2))*exp(-(x-v1-v2*t1)^2/(2*exp(2*v3+2*v4*t2)))}
compare("d",dnorm(1,4+5*2,exp(6+7*3)),f(1,2,3,4,5,6,7))
norm_p12_fd=Deriv(f,c("v1","v2","v3","v4"),nderiv=1)
norm_p12_fdd=Deriv(f,c("v1","v2","v3","v4"),nderiv=2)
#
p=function(x,t1,t2,v1,v2,v3,v4){pnorm((x-v1-v2*t1)/exp(v3+v4*t2))}
compare("p",pnorm(1,4+5*2,exp(6+7*3)),p(1,2,3,4,5,6,7))
norm_p12_pd=Deriv(p,c("v1","v2","v3","v4"),nderiv=1)
norm_p12_pdd=Deriv(p,c("v1","v2","v3","v4"),nderiv=2)
#
logf=function(x,t1,t2,v1,v2,v3,v4){-0.5*log(2*pi)-v3-v4*t2-(x-v1-v2*t1)^2/(2*exp(2*v3+2*v4*t2))}
compare("l",dnorm(1,4+5*2,exp(6+7*3),log=TRUE),logf(1,2,3,4,5,6,7))
norm_p12_logfdd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=2)
norm_p12_logfddd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=3)
#
sink("80c_norm_p12_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_fd=")
print.function(norm_p12_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_fdd=")
print.function(norm_p12_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_pd=")
print.function(norm_p12_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_pdd=")
print.function(norm_p12_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_logfdd=")
print.function(norm_p12_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p12_logfddd=")
print.function(norm_p12_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_f1fa=function(x,t1,t2,v1,v2,v3,v4){
	vf=Vectorize(norm_p12_fd)
	f1=vf(x,t1,t2,v1,v2,v3,v4)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_f2fa=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_fdd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_p1fa=function(x,t1,t2,v1,v2,v3,v4){
	vf=Vectorize(norm_p12_pd)
	p1=vf(x,t1,t2,v1,v2,v3,v4)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_p2fa=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_pdd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	p2=deriv_copyfdd(temp1,nx,dim=4)
	return(p2)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_mu1fa=function(alpha,t1,t2,v1,v2,v3,v4){
	x=qnorm((1-alpha),mean=v1+v2*t1,sd=exp(v3+v4*t2))
	vf=Vectorize(norm_p12_pd)
	mu1=-vf(x,t1,t2,v1,v2,v3,v4)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_mu2fa=function(alpha,t1,t2,v1,v2,v3,v4){
	x=qnorm((1-alpha),mean=v1+v2*t1,sd=exp(v3+v4*t2))
	nx=length(x)
	vf=Vectorize(norm_p12_pdd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	mu2=-deriv_copyfdd(temp1,nx,dim=4)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_ldda=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_logfdd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p12_lddda=function(x,t1,t2,v1,v2,v3,v4){
	nx=length(x)
	vf=Vectorize(norm_p12_logfddd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
