#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# note that in this function definition, v1=lambda, v2=mu, v3=sigma
# which is the ordering in R
# which is different from my main code...be careful
# the switch occurs below when I call the Deriv functions
f=function(x,v1,v2,v3){log(v1)-log(v3)-(1+v1)*log((x-v2)/v3)-((x-v2)/v3)^(-v1)}
frechet_k1_logfdd=Deriv(f,c("v1","v3"),nderiv=2)
frechet_k1_logfddd=Deriv(f,c("v1","v3"),nderiv=3)
#
sink("51c_frechet_k1_derivs.R")
#
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
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"frechet_k1_ldda=function(x,v1,v2,kloc){
# the v1 coming in here is sigma, and the v2 is lambda, following my pu code
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(frechet_k1_logfdd)
	temp1=vf(x,v2,kloc,v1) #these are in lambda, mu, sigma order
	temp2=apply(temp1,1,sum)/nx
	for (i in 1:2){
		for (j in 1:2){
			ldd[3-i,3-j]=temp2[(i-1)*2+j] #I have to switch the params around
		}
	}
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
# I have to switch because my pu code orders sigma and lambda differently
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(frechet_k1_logfddd)
	temp1=vf(x,v2,kloc,v1) #these are in lambda, mu, sigma order
	temp2=apply(temp1,1,sum)/nx
	for (i in 1:2){
		for (j in 1:2){
			for (k in 1:2){
				lddd[3-i,3-j,3-k]=temp2[(i-1)*4+(j-1)*2+k] #I have to switch the params around
			}
		}
	}
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
