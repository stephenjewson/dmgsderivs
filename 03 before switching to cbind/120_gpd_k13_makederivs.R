#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# note that in this function definition, v1=mu, v2=sigma, v3=xi
# which is the ordering in R
# my code uses v1=sigma, v2=xi
f=function(x,v1,v2,v3){(1/v2)*(1+v3*((x-v1)/v2))^(-((v3+1)/v3))}
compare("d",dgpd(5,2,3,4),f(5,2,3,4))
gpd_k13_fd=Deriv(f,c("v2"),nderiv=1)
gpd_k13_fdd=Deriv(f,c("v2"),nderiv=2)
#
p=function(x,v1,v2,v3){1-(1+v3*((x-v1)/v2))^(-1/v3)}
compare("p",pgpd(5,2,3,4),p(5,2,3,4))
gpd_k13_pd=Deriv(p,c("v2"),nderiv=1)
gpd_k13_pdd=Deriv(p,c("v2"),nderiv=2)
#
logf=function(x,v1,v2,v3){-log(v2)-((1+v3)/v3)*log(1+v3*(x-v1)/v2)}
compare("l",dgpd(5,2,3,4,log=TRUE),logf(5,2,3,4))
gpd_k13_logfdd=Deriv(logf,c("v2"),nderiv=2)
gpd_k13_logfddd=Deriv(logf,c("v2"),nderiv=3)
#
sink("120c_gpd_k13_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k13_fd=")
print.function(gpd_k13_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k13_fdd=")
print.function(gpd_k13_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k13_pd=")
print.function(gpd_k13_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k13_pdd=")
print.function(gpd_k13_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k13_logfdd=")
print.function(gpd_k13_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k13_logfddd=")
print.function(gpd_k13_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k13_f1fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code

	nx=length(x)
	f1=matrix(0,1,nx)
	vf=Vectorize(gpd_k13_fd)
	v2=movexiawayfromzero(v2)
	f1[1,]=vf(x,kloc,v1,v2)

	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k13_f2fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code

	nx=length(x)
	f2=array(0,c(1,1,nx))
	vf=Vectorize(gpd_k13_fdd)
	v2=movexiawayfromzero(v2)
	f2[1,1,]=vf(x,kloc,v1,v2)

	return(f2)
}\n"
)
###cat("############################################################\n")
####
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gpd_k13_p1fa=function(x,v1,v2,kloc){
#### the v1 coming in here is sigma, and the v2 is lambda, following my pu code
###
###	nx=length(x)
###	p1=matrix(0,1,nx)
###	vf=Vectorize(gpd_k13_pd)
###	v2=movexiawayfromzero(v2)
###	p1[1,]=vf(x,kloc,v1,v2)
###
###	return(p1)
###}\n"
###)
###cat("############################################################\n")
####
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"gpd_k13_p2fa=function(x,v1,v2,kloc){
#### the v1 coming in here is sigma, and the v2 is lambda, following my pu code
###
###	nx=length(x)
###	p2=array(0,c(1,1,nx))
###	vf=Vectorize(gpd_k13_pdd)
###	v2=movexiawayfromzero(v2)
###	p2[1,1,]=vf(x,kloc,v1,v2)
###
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k13_mu1fa=function(alpha,v1,v2,kloc){
	x=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
	nx=length(x)
	mu1=array(0,c(1,nx))
	vf=Vectorize(gpd_k13_pd)
	mu1[1,]=-vf(x,v1,v2,kloc)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k13_mu2fa=function(alpha,v1,v2,kloc){
	x=extraDistr::qgpd((1-alpha),mu=kloc,sigma=v1,xi=v2)
	nx=length(x)
	mu2=array(0,c(1,1,nx))
	vf=Vectorize(gpd_k13_pdd)
	mu2[1,1,]=-vf(x,v1,v2,kloc)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k13_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code

	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(gpd_k13_logfdd)
	v2=movexiawayfromzero(v2)
	ldd[1,1]=sum(vf(x,kloc,v1,v2))/nx

	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k13_lddda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch because my pu code orders sigma and lambda differently

	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(gpd_k13_logfddd)
	v2=movexiawayfromzero(v2)
	lddd[1,1,1]=sum(vf(x,kloc,v1,v2))/nx

	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
