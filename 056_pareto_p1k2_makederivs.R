#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
library(actuar)
#
# v1=shape (a)
# v2=scale (b)
# note that the actuar version of the pareto is different
#f=function(x,t,v1,v2,v3){v1*(v2^v1)/((x)^(v1+1))}
f=function(x,t,v1,v2,v3){exp(-v1-v2*t)*(v3^exp(-v1-v2*t))/((x)^(exp(-v1-v2*t)+1))}
compare("d",extraDistr::dpareto(5,exp(-3-2*4),2),f(5,4,3,2,2))
pareto_p1k2_fd=Deriv(f,c("v1","v2"),nderiv=1)
pareto_p1k2_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
#p=function(x,t,v1,v2,v3){1-(v2/x)^v1}
p=function(x,t,v1,v2,v3){1-(v3/x)^exp(-v1-v2*t)}
compare("p",extraDistr::ppareto(5,exp(-3-2*4),2),p(5,4,3,2,2))
pareto_p1k2_pd=Deriv(p,c("v1","v2"),nderiv=1)
pareto_p1k2_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
#logf=function(x,t,v1,v2,v3){log(v1)+v1*log(v2)-(v1+1)*log(x)}
logf=function(x,t,v1,v2,v3){-v1-v2*t+exp(-v1-v2*t)*log(v3)-(exp(-v1-v2*t)+1)*log(x)}
compare("l",extraDistr::dpareto(5,exp(-3-2*4),2,log=TRUE),logf(5,4,3,2,2))
pareto_p1k2_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
pareto_p1k2_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("056c_pareto_p1k2_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k2_fd=")
print.function(pareto_p1k2_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k2_fdd=")
print.function(pareto_p1k2_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k2_pd=")
print.function(pareto_p1k2_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k2_pdd=")
print.function(pareto_p1k2_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k2_logfdd=")
print.function(pareto_p1k2_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k2_logfddd=")
print.function(pareto_p1k2_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_p1k2_f1fa=function(x,t,v1,v2,kscale){
	vf=Vectorize(pareto_p1k2_fd,\"x\")
	f1=vf(x,t,v1,v2,kscale)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_p1k2_f2fa=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_fdd,\"x\")
	temp1=vf(x,t,v1,v2,kscale)
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
"pareto_p1k2_p1fa=function(x,t,v1,v2,kscale){
	vf=Vectorize(pareto_p1k2_pd,\"x\")
	p1=vf(x,t,v1,v2,kscale)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_p1k2_p2fa=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_pdd,\"x\")
	temp1=vf(x,t,v1,v2,kscale)
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
"pareto_p1k2_mu1fa=function(alpha,t,v1,v2,kscale){
	x=extraDistr::qpareto((1-alpha),a=exp(-v1-v2*t),b=kscale)
	vf=Vectorize(pareto_p1k2_pd,\"x\")
	mu1=-vf(x,t,v1,v2,kscale)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_p1k2_mu2fa=function(alpha,t,v1,v2,kscale){
	x=extraDistr::qpareto((1-alpha),a=exp(-v1-v2*t),b=kscale)
	nalpha=length(alpha)
	vf=Vectorize(pareto_p1k2_pdd,\"x\")
	temp1=vf(x,t,v1,v2,kscale)
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
"pareto_p1k2_ldda=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_logfdd,\"x\")
	temp1=vf(x,t,v1,v2,kscale)
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
"pareto_p1k2_lddda=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k2_logfddd,\"x\")
	temp1=vf(x,t,v1,v2,kscale)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))

