#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2,v3){(1/sqrt(2*pi))*(1/v3)*exp(-(x-v1-v2*t)^2/(2*v3*v3))}
compare("d",dnorm(1,3+4*2,5),f(1,2,3,4,5))
norm_p1_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
norm_p1_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
p=function(x,t,v1,v2,v3){pnorm((x-v1-v2*t)/v3)}
compare("p",pnorm(1,3+4*2,5),p(1,2,3,4,5))
norm_p1_pd=Deriv(p,c("v1","v2","v3"),nderiv=1)
norm_p1_pdd=Deriv(p,c("v1","v2","v3"),nderiv=2)
#
logf=function(x,t,v1,v2,v3){-0.5*log(2*pi)-log(v3)-(x-v1-v2*t)^2/(2*v3*v3)}
compare("l",dnorm(1,3+4*2,5,log=TRUE),logf(1,2,3,4,5))
norm_p1_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
norm_p1_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("60c_norm_p1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p1_fd=")
print.function(norm_p1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p1_fdd=")
print.function(norm_p1_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p1_pd=")
print.function(norm_p1_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p1_pdd=")
print.function(norm_p1_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p1_logfdd=")
print.function(norm_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("norm_p1_logfddd=")
print.function(norm_p1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_f1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(norm_p1_fd)
	f1=vf(x,t,v1,v2,v3)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_f2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_fdd)
	temp1=vf(x,t,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_p1fa=function(x,t,v1,v2,v3){
	vf=Vectorize(norm_p1_pd)
	p1=vf(x,t,v1,v2,v3)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_p2fa=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_pdd)
	temp1=vf(x,t,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_mu1fa=function(alpha,t,v1,v2,v3){
	x=qnorm((1-alpha),mean=v1+v2*t,sd=v3)
	vf=Vectorize(norm_p1_pd)
	mu1=-vf(x,t,v1,v2,v3)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_mu2fa=function(alpha,t,v1,v2,v3){
	x=qnorm((1-alpha),mean=v1+v2*t,sd=v3)
	nx=length(x)
	vf=Vectorize(norm_p1_pdd)
	temp1=vf(x,t,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_ldda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_logfdd)
	temp1=vf(x,t,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"norm_p1_lddda=function(x,t,v1,v2,v3){
	nx=length(x)
	vf=Vectorize(norm_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
