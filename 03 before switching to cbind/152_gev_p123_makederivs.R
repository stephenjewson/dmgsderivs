#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t1,t2,t3,v1,v2,v3,v4,v5,v6){
	(exp(-v3-v4*t2))*((1+(v5+v6*t3)*((x-v1-v2*t1)*exp(-v3-v4*t2)))^(-1/(v5+v6*t3)-1))*
		exp(-((1+(v5+v6*t3)*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/(v5+v6*t3))))
}
compare("d",extraDistr::dgev(1,5+.6*2,exp(.7+.8*3),0.9+0.09*4),f(1,2,3,4,5,.6,.7,.8,.9,.09))
gev_p123_fd=Deriv(f,c("v1","v2","v3","v4","v5","v6"),nderiv=1)
gev_p123_fdd=Deriv(f,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
#
p=function(x,t1,t2,t3,v1,v2,v3,v4,v5,v6){
	exp(-(1+(v5+v6*t3)*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/(v5+v6*t3)))
	}
compare("p",extraDistr::pgev(1,5+.6*2,exp(.7+.8*3),0.9+0.09*4),p(1,2,3,4,5,.6,.7,.8,.9,.09))
gev_p123_pd=Deriv(p,c("v1","v2","v3","v4","v5","v6"),nderiv=1)
gev_p123_pdd=Deriv(p,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
#
logf=function(x,t1,t2,t3,v1,v2,v3,v4,v5,v6){
	-v3-v4*t2-(1+1/(v5+v6*t3))*log(1+(v5+v6*t3)*(x-v1-v2*t1)*exp(-v3-v4*t2))-
		(1+(v5+v6*t3)*(x-v1-v2*t1)*exp(-v3-v4*t2))^(-1/(v5+v6*t3))
}
compare("l",extraDistr::dgev(1,5+.6*2,exp(.7+.8*3),0.9+0.09*4,log=TRUE),logf(1,2,3,4,5,.6,.7,.8,.9,.09))
gev_p123_logfdd=Deriv(logf,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
gev_p123_logfddd=Deriv(logf,c("v1","v2","v3","v4","v5","v6"),nderiv=3)
#
sink("152c_gev_p123_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_fd=")
print.function(gev_p123_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_fdd=")
print.function(gev_p123_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_pd=")
print.function(gev_p123_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_pdd=")
print.function(gev_p123_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_logfdd=")
print.function(gev_p123_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_logfddd=")
print.function(gev_p123_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_f1fa=function(x,t01,t02,t03,v1,v2,v3,v4,v5,v6){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p123_fd)
	f1=vf(x,t01,t02,t03,v1,v2,v3,v4,v5,v6)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_f2fa=function(x,t01,t02,t03,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p123_fdd)
	temp1=vf(x,t01,t02,t03,v1,v2,v3,v4,v5,v6)
	f2=deriv_copyfdd(temp1,nx,dim=6)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p123_p1fa=function(x,t01,t02,t03,v1,v2,v3,v4,v5,v6){
###
###	v3=movexiawayfromzero(v3)
###
###	vf=Vectorize(gev_p123_pd)
###	p1=vf(x,t01,t02,t03,v1,v2,v3,v4,v5,v6)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p123_p2fa=function(x,t01,t02,t03,v1,v2,v3,v4,v5,v6){
###	nx=length(x)
###
###	v3=movexiawayfromzero(v3)
###
###	vf=Vectorize(gev_p123_pdd)
###	temp1=vf(x,t01,t02,t03,v1,v2,v3,v4,v5,v6)
###	p2=deriv_copyfdd(temp1,nx,dim=6)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_mu1fa=function(alpha,t01,t02,t03,v1,v2,v3,v4,v5,v6){
	x=qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=v5+v6*t03)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p123_pd)
	mu1=-vf(x,t01,t02,t03,v1,v2,v3,v4,v5,v6)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_mu2fa=function(alpha,t01,t02,t03,v1,v2,v3,v4,v5,v6){
	x=qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=v5+v6*t03)
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p123_pdd)
	temp1=vf(x,t01,t02,t03,v1,v2,v3,v4,v5,v6)
	mu2=-deriv_copyfdd(temp1,nx,dim=6)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_ldda=function(x,t1,t2,t3,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p123_logfdd)
	temp1=vf(x,t1,t2,t3,v1,v2,v3,v4,v5,v6)
	ldd=deriv_copyldd(temp1,nx,dim=6)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_lddda=function(x,t1,t2,t3,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p123_logfddd)
	temp1=vf(x,t1,t2,t3,v1,v2,v3,v4,v5,v6)
	lddd=deriv_copylddd(temp1,nx,dim=6)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
