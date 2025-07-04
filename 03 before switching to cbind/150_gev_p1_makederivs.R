#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2,v3,v4){
	(1/v3)*((1+v4*((x-v1-v2*t)/v3))^(-1/v4-1))*exp(-((1+v4*(x-v1-v2*t)/v3)^(-1/v4)))
}
compare("d",extraDistr::dgev(1,3+4*2,5,0.1),f(1,2,3,4,5,0.1))
gev_p1_fd=Deriv(f,c("v1","v2","v3","v4"),nderiv=1)
gev_p1_fdd=Deriv(f,c("v1","v2","v3","v4"),nderiv=2)
#
p=function(x,t,v1,v2,v3,v4){
	exp(-((1+v4*(x-v1-v2*t)/v3)^(-1/v4)))
	}
compare("p",extraDistr::pgev(1,3+4*2,5,0.1),p(1,2,3,4,5,0.1))
gev_p1_pd=Deriv(p,c("v1","v2","v3","v4"),nderiv=1)
gev_p1_pdd=Deriv(p,c("v1","v2","v3","v4"),nderiv=2)
#
logf=function(x,t,v1,v2,v3,v4){
	-log(v3)-(1+1/v4)*log(1+v4*(x-v1-v2*t)/v3)-(1+v4*(x-v1-v2*t)/v3)^(-1/v4)
}
compare("l",extraDistr::dgev(1,3+4*2,5,0.1,log=TRUE),logf(1,2,3,4,5,0.1))
gev_p1_logfdd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=2)
gev_p1_logfddd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=3)
#
sink("150c_gev_p1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p1_fd=")
print.function(gev_p1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p1_fdd=")
print.function(gev_p1_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p1_pd=")
print.function(gev_p1_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p1_pdd=")
print.function(gev_p1_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p1_logfdd=")
print.function(gev_p1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p1_logfddd=")
print.function(gev_p1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_f1fa=function(x,t,v1,v2,v3,v4){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_fd)
	f1=vf(x,t,v1,v2,v3,v4)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_f2fa=function(x,t,v1,v2,v3,v4){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_fdd)
	temp1=vf(x,t,v1,v2,v3,v4)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1_p1fa=function(x,t,v1,v2,v3,v4){
###
###	v3=movexiawayfromzero(v3)
###
###	vf=Vectorize(gev_p1_pd)
###	p1=vf(x,t,v1,v2,v3,v4)
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gev_p1_p2fa=function(x,t,v1,v2,v3,v4){
###	nx=length(x)
###
###	v3=movexiawayfromzero(v3)
###
###	vf=Vectorize(gev_p1_pdd)
###	temp1=vf(x,t,v1,v2,v3,v4)
###	p2=deriv_copyfdd(temp1,nx,dim=4)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_mu1fa=function(alpha,t,v1,v2,v3,v4){
	x=qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=v4)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_pd)
	mu1=-vf(x,t,v1,v2,v3,v4)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_mu2fa=function(alpha,t,v1,v2,v3,v4){
	x=qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=v4)
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_pdd)
	temp1=vf(x,t,v1,v2,v3,v4)
	mu2=-deriv_copyfdd(temp1,nx,dim=4)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_ldda=function(x,t,v1,v2,v3,v4){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_logfdd)
	temp1=vf(x,t,v1,v2,v3,v4)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_lddda=function(x,t,v1,v2,v3,v4){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3,v4)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
