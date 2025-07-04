#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# note that in this function definition, v1=mu, v2,3=sigma, v4=lambda
# which is different from my main code...be careful
f=function(x,t,v1,v2,v3,v4){log(v4)-v2-t*v3-(1+v4)*log((x-v1)/exp(v2+t*v3))-((x-v1)/exp(v2+t*v3))^(-v4)}
frechet_p2k1_logfdd=Deriv(f,c("v2","v3","v4"),nderiv=2)
frechet_p2k1_logfddd=Deriv(f,c("v2","v3","v4"),nderiv=3)
#
sink("71c_frechet_p2k1_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_logfdd=")
print.function(frechet_p2k1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("frechet_p2k1_logfddd=")
print.function(frechet_p2k1_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_ldda=function(x,t,v1,v2,v3,kloc){
# the v1,2 coming in here is sigma, and the v3 is lambda, following my pu code
# I have to switch below
	nx=length(x)
	ldd=matrix(0,3,3)
	vf=Vectorize(frechet_p2k1_logfdd)
	temp1=vf(x,t,kloc,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_p2k1_lddda=function(x,t,v1,v2,v3,kloc){
# the v1,2 coming in here is sigma, and the v3 is lambda, following my pu code
# I have to switch below
	nx=length(x)
	lddd=array(0,c(3,3,3))
	vf=Vectorize(frechet_p2k1_logfddd)
	temp1=vf(x,t,kloc,v1,v2,v3) #these are in mu, sigma, lambda order
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
