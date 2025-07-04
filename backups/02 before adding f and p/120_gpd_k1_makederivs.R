#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# note that in this function definition, v1=mu, v2=sigma, v3=xi
# which is the ordering in R
# my code uses v1=sigma, v2=xi
f=function(x,v1,v2,v3){-log(v2)-((1+v3)/v3)*log(1+v3*(x-v1)/v2)}
gpd_k1_logfdd=Deriv(f,c("v2","v3"),nderiv=2)
gpd_k1_logfddd=Deriv(f,c("v2","v3"),nderiv=3)
#
sink("120c_gpd_k1_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k1_logfdd=")
print.function(gpd_k1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gpd_k1_logfddd=")
print.function(gpd_k1_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k1_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(gpd_k1_logfdd)

	v2=movexiawayfromzero(v2)

	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gpd_k1_lddda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
# I have to switch because my pu code orders sigma and lambda differently
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(gpd_k1_logfddd)

	v2=movexiawayfromzero(v2)

	temp1=vf(x,kloc,v1,v2) #these are in mu, sigma, xi order
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
