#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,t,v1,v2,v3,v4){-log(v3)-(1+1/v4)*log(1+v4*(x-v1-t*v2)/v3)-(1+v4*(x-v1-t*v2)/v3)^(-1/v4)}
gev_p1_logfdd=Deriv(f,c("v1","v2","v3","v4"),nderiv=2)
gev_p1_logfddd=Deriv(f,c("v1","v2","v3","v4"),nderiv=3)
#
sink("150c_gev_p1_derivs.R")
#
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
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p1_ldda=function(x,t,v1,v2,v3,v4){
	nx=length(x)
	ldd=matrix(0,4,4)

	v4=movexiawayfromzero(v4)

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
	lddd=array(0,c(4,4,4))

	v4=movexiawayfromzero(v4)

	vf=Vectorize(gev_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3,v4)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
