#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# note that in this function definition, v1=mu, v2=sigma, v3=lambda
# which is different from my main code...be careful
f=function(x,v1,v2,v3){(v3/v2)*((x-v1)/v2)^(-1-v3)*exp(-((x-v1)/v2)^(-v3))}
compare("d",dfrechet(5,2,3,4),f(5,3,4,2))
frechet_k1_fd=Deriv(f,c("v2","v3"),nderiv=1)
frechet_k1_fdd=Deriv(f,c("v2","v3"),nderiv=2)
#
p=function(x,v1,v2,v3){exp(-((x-v1)/v2)^(-v3))}
compare("p",pfrechet(5,2,3,4),p(5,3,4,2))
frechet_k1_pd=Deriv(p,c("v2","v3"),nderiv=1)
frechet_k1_pdd=Deriv(p,c("v2","v3"),nderiv=2)
#
logf=function(x,v1,v2,v3){log(v3)-log(v2)-(1+v3)*log((x-v1)/v2)-((x-v1)/v2)^(-v3)}
compare("l",dfrechet(5,2,3,4,log=TRUE),logf(5,3,4,2))
frechet_k1_logfdd=Deriv(logf,c("v2","v3"),nderiv=2)
frechet_k1_logfddd=Deriv(logf,c("v2","v3"),nderiv=3)
#
sink("51c_frechet_k1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_k1_fd=")
print.function(frechet_k1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_k1_fdd=")
print.function(frechet_k1_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_k1_pd=")
print.function(frechet_k1_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_k1_pdd=")
print.function(frechet_k1_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_k1_logfdd=")
print.function(frechet_k1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_k1_logfddd=")
print.function(frechet_k1_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_k1_f1fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch below
	nx=length(x)
	f1=matrix(0,2,nx)
	vf=Vectorize(frechet_k1_fd)
	temp1=vf(x,kloc,v1,v2)
	f1=deriv_copyldd(temp1,nx,dim=2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_k1_f2fa=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch below
	nx=length(x)
	f2=array(0,c(2,2,nx))
	vf=Vectorize(frechet_k1_fdd)
	temp1=vf(x,kloc,v1,v2)
	f2=deriv_copyldd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_k1_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch below
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(frechet_k1_logfdd)
	temp1=vf(x,kloc,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_k1_lddda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch below
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(frechet_k1_logfddd)
	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, lambda order
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
