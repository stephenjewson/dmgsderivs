#
# make derivative codes for fitdistpu, one model at a time
#
# for exp, I don't actually use the mu and p routines
# but I add them here, so that I can derive them for exp_p1, where I do use them
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1){v1*exp(-v1*x)}
compare("d",dexp(1,2),f(1,2))
exp_fd=Deriv(f,"v1",nderiv=1)
exp_fdd=Deriv(f,"v1",nderiv=2)
#
p=function(x,v1){1-exp(-v1*x)}
compare("p",pexp(1,2),p(1,2))
exp_pd=Deriv(p,"v1",nderiv=1)
exp_pdd=Deriv(p,"v1",nderiv=2)
#
logf=function(x,v1){log(v1)-v1*x}
compare("l",dexp(1,2,log=TRUE),logf(1,2))
exp_logfdd=Deriv(logf,"v1",nderiv=2)
exp_logfddd=Deriv(logf,"v1",nderiv=3)
#
sink("10c_exp_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_fd=")
print.function(exp_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_fdd=")
print.function(exp_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_pd=")
print.function(exp_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_pdd=")
print.function(exp_pdd)
cat("######################################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_logfdd=")
print.function(exp_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_logfddd=")
print.function(exp_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"exp_f1fa=function(x,v1){
	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(exp_fd)
	f1[1,]=vf(x,v1)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"exp_f2fa=function(x,v1){
	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(exp_fdd)
	f2[1,1,]=vf(x,v1)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p1fa=function(x,v1){
	nx=length(x)
	p1=matrix(0,1,nx)
	vf=Vectorize(exp_pd)
	p1[1,]=vf(x,v1)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @inheritParams manf\n")
cat(
"exp_p2fa=function(x,v1){
	nx=length(x)
	p2=array(0,c(1,1,nx))
	vf=Vectorize(exp_pdd)
	p2[1,1,]=vf(x,v1)
	return(p2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' Minus the first derivative of the cdf, at alpha\n")
###cat("#' @inheritParams manf\n")
###cat(
###"exp_mu1fa=function(alpha,v1){
###	x=qexp((1-alpha),rate=v1)
###	vf=Vectorize(exp_pd)
###	mu1=-vf(x,v1)
###	return(mu1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' Minus the second derivative of the cdf, at alpha\n")
###cat("#' @inheritParams manf\n")
###cat(
###"exp_mu2fa=function(alpha,v1){
###	x=qexp((1-alpha),rate=v1)
###	nx=length(x)
###	vf=Vectorize(exp_pdd)
###	mu2=-vf(x,v1)
###	return(mu2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"exp_ldda=function(x,v1){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(exp_logfdd)
	ldd[1,1]=sum(vf(x,v1))/nx
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"exp_lddda=function(x,v1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(exp_logfddd)
	lddd[1,1,1]=sum(vf(x,v1))/nx
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

