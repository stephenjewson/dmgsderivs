#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
#
# note that in this function definition, v1=mu, v2,3=sigma, v4=lambda
# which is different from my main code...be careful
f=function(x,t,v1,v2,v3,v4){(v4/exp(v2+v3*t))*((x-v1)/exp(v2+v3*t))^(-1-v4)*exp(-((x-v1)/exp(v2+v3*t))^(-v4))}
compare("d",dfrechet(10,2,0,exp(.1+2*.3)),f(10,.3,0,.1,2,2))
frechet_p2k1_fd=Deriv(f,c("v2","v3","v4"),nderiv=1)
frechet_p2k1_fdd=Deriv(f,c("v2","v3","v4"),nderiv=2)
#
p=function(x,t,v1,v2,v3,v4){exp(-((x-v1)/exp(v2+v3*t))^(-v4))}
compare("p",pfrechet(10,2,0,exp(.1+2*.3)),p(10,.3,0,.1,2,2))
frechet_p2k1_pd=Deriv(p,c("v2","v3","v4"),nderiv=1)
frechet_p2k1_pdd=Deriv(p,c("v2","v3","v4"),nderiv=2)
#
logf=function(x,t,v1,v2,v3,v4){log(v4)-v2-v3*t-(1+v4)*log((x-v1))+(1+v4)*(v2+v3*t)-((x-v1)*exp(-v2-v3*t))^(-v4)}
compare("l",dfrechet(10,2,0,exp(.1+2*.3),log=TRUE),logf(10,.3,0,.1,2,2))
frechet_p2k1_logfdd=Deriv(logf,c("v2","v3","v4"),nderiv=2)
frechet_p2k1_logfddd=Deriv(logf,c("v2","v3","v4"),nderiv=3)
#
sink("071c_frechet_p2k1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_fd=")
print.function(frechet_p2k1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_fdd=")
print.function(frechet_p2k1_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_pd=")
print.function(frechet_p2k1_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_pdd=")
print.function(frechet_p2k1_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_logfdd=")
print.function(frechet_p2k1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_logfddd=")
print.function(frechet_p2k1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_f1fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_p2k1_fd,\"x\")
	f1=vf(x,t0,kloc,v1,v2,v3)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_f2fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_fdd,\"x\")
	temp1=vf(x,t0,kloc,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The first derivative of the cdf\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_p1fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	vf=Vectorize(frechet_p2k1_pd,\"x\")
	p1=vf(x,t0,kloc,v1,v2,v3)
	return(p1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the cdf\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_p2fa=function(x,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_pdd,\"x\")
	temp1=vf(x,t0,kloc,v1,v2,v3)
	p2=deriv_copyfdd(temp1,nx,dim=3)
	return(p2)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the first derivative of the cdf, at alpha\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_mu1fa=function(alpha,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	x=qfrechet((1-alpha),mu=kloc,sigma=exp(v1+v2*t0),lambda=v3)
#	x=qfrechet((1-alpha),mu=kloc,sigma=exp(v2+v3*t0),lambda=vkloc)
	vf=Vectorize(frechet_p2k1_pd,\"x\")
	mu1=-vf(x,t0,kloc,v1,v2,v3)
	return(mu1)
}\n"
)
cat("############################################################\n")
#
cat("#' Minus the second derivative of the cdf, at alpha\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_mu2fa=function(alpha,t0,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	x=qfrechet((1-alpha),mu=kloc,sigma=exp(v1+v2*t0),lambda=v3)
	nx=length(x)
	vf=Vectorize(frechet_p2k1_pdd,\"x\")
	temp1=vf(x,t0,kloc,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_ldda=function(x,t,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_logfdd,c(\"x\",\"t\"))
	temp1=vf(x,t,kloc,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_lddda=function(x,t,v1,v2,v3,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my cp code
# I have to switch below
	nx=length(x)
	vf=Vectorize(frechet_p2k1_logfddd,c(\"x\",\"t\"))
	temp1=vf(x,t,kloc,v1,v2,v3) #these are in mu, sigma, lambda order
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
