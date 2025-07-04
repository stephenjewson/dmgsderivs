#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(extraDistr)
library(actuar)
#
# v1=shape (a)
# v2=scale (b)
# note that the actuar version of the pareto is different
f=function(x,v1,v2){v1*(v2^v1)/((x)^(v1+1))}
compare("d",extraDistr::dpareto(4,3,2),f(4,3,2))
pareto_k2_fd=Deriv(f,"v1",nderiv=1)
pareto_k2_fdd=Deriv(f,"v1",nderiv=2)
#
p=function(x,v1,v2){1-(v2/x)^v1}
compare("p",extraDistr::ppareto(4,3,2),p(4,3,2))
pareto_k2_pd=Deriv(p,"v1",nderiv=1)
pareto_k2_pdd=Deriv(p,"v1",nderiv=2)
#
logf=function(x,v1,v2){log(v1)+v1*log(v2)-(v1+1)*log(x)}
compare("l",extraDistr::dpareto(4,3,2,log=TRUE),logf(4,3,2))
pareto_k2_logfdd=Deriv(logf,"v1",nderiv=2)
pareto_k2_logfddd=Deriv(logf,"v1",nderiv=3)
#
sink("11c_pareto_k2_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k2_fd=")
print.function(pareto_k2_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k2_fdd=")
print.function(pareto_k2_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k2_pd=")
print.function(pareto_k2_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k2_pdd=")
print.function(pareto_k2_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k2_logfdd=")
print.function(pareto_k2_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k2_logfddd=")
print.function(pareto_k2_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_f1fa=function(x,v1,kscale){
	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(pareto_k2_fd)
	f1[1,]=vf(x,v1,v2=kscale)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_f2fa=function(x,v1,kscale){
	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(pareto_k2_fdd)
	f2[1,1,]=vf(x,v1,v2=kscale)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_p1fa=function(x,v1,kscale){
	nx=length(x)
	p1=matrix(0,1,nx)
	vf=Vectorize(pareto_k2_pd)
	p1[1,]=vf(x,v1,v2=kscale)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_p2fa=function(x,v1,kscale){
	nx=length(x)
	p2=array(0,c(1,1,nx))
	vf=Vectorize(pareto_k2_pdd)
	p2[1,1,]=vf(x,v1,v2=kscale)
	return(p2)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_mu1fa=function(alpha,v1,kscale){
	x=extraDistr::qpareto((1-alpha),a=v1,b=kscale)
	vf=Vectorize(pareto_k2_pd)
	mu1=-vf(x,v1,kscale)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_mu2fa=function(alpha,v1,kscale){
	x=qpareto((1-alpha),a=v1,b=kscale)
	nx=length(x)
	vf=Vectorize(pareto_k2_pdd)
	mu2=-vf(x,v1,kscale)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_ldda=function(x,v1,kscale){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(pareto_k2_logfdd)
	ldd[1,1]=sum(vf(x,v1,v2=kscale))/nx
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k2_lddda=function(x,v1,kscale){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	temp=0
	vf=Vectorize(pareto_k2_logfddd)
	lddd[1,1,1]=sum(vf(x,v1,v2=kscale))/nx
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

