#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1,v2,v3){(1/v2)*((1+v3*((x-v1)/v2))^(-1/v3-1))*exp(-((1+v3*(x-v1)/v2)^(-1/v3)))}
compare("d",extraDistr::dgev(1,2,3,0.1),f(1,2,3,0.1))
gev_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
gev_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
logf=function(x,v1,v2,v3){-log(v2)-(1+1/v3)*log(1+v3*(x-v1)/v2)-(1+v3*(x-v1)/v2)^(-1/v3)}
compare("l",extraDistr::dgev(1,2,3,0.1,log=TRUE),logf(1,2,3,0.1))
gev_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
gev_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("110c_gev_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_fd=")
print.function(gev_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_fdd=")
print.function(gev_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_logfdd=")
print.function(gev_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_logfddd=")
print.function(gev_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_ldda=function(x,v1,v2,v3){
	nx=length(x)
	ldd=matrix(0,3,3)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfdd)
	temp1=vf(x,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_lddda=function(x,v1,v2,v3){
	nx=length(x)
	lddd=array(0,c(3,3,3))

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfddd)
	temp1=vf(x,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
