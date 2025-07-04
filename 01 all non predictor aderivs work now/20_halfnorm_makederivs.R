#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1){log(2)+log(v1)-log(pi)-x*x*v1*v1/pi}
halfnorm_logfdd=Deriv(f,"v1",nderiv=2)
halfnorm_logfddd=Deriv(f,"v1",nderiv=3)
#
sink("20c_halfnorm_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("halfnorm_logfdd=")
print.function(halfnorm_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("halfnorm_logfddd=")
print.function(halfnorm_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"halfnorm_ldda=function(x,v1){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(halfnorm_logfdd)
	ldd[1,1]=sum(vf(x,v1))/nx
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"halfnorm_lddda=function(x,v1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(halfnorm_logfddd)
	lddd[1,1,1]=sum(vf(x,v1))/nx
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

