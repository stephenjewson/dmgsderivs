#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(extraDistr)
#
f=function(x,v1,v2,v3){(1/v2)*((1+v3*((x-v1)/v2))^(-1/v3-1))*exp(-((1+v3*(x-v1)/v2)^(-1/v3)))}
compare("d",extraDistr::dgev(1,2,3,0.1),f(1,2,3,0.1))
gev_k3_fd=Deriv(f,c("v1","v2"),nderiv=1)
gev_k3_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,v1,v2,v3){exp(-((1+v3*(x-v1)/v2)^(-1/v3)))}
compare("p",extraDistr::pgev(1,2,3,0.1),p(1,2,3,0.1))
gev_k3_pd=Deriv(p,c("v1","v2"),nderiv=1)
gev_k3_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2,v3){-log(v2)-(1+1/v3)*log(1+v3*(x-v1)/v2)-(1+v3*(x-v1)/v2)^(-1/v3)}
compare("l",extraDistr::dgev(1,2,3,0.1,log=TRUE),logf(1,2,3,0.1))
gev_k3_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
gev_k3_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("53c_gev_k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_k3_fd=")
print.function(gev_k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_k3_fdd=")
print.function(gev_k3_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_k3_pd=")
print.function(gev_k3_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_k3_pdd=")
print.function(gev_k3_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_k3_logfdd=")
print.function(gev_k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_k3_logfddd=")
print.function(gev_k3_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gev_k3_f1fa=function(x,v1,v2,kshape){
	nx=length(x)
	f1=matrix(0,2,nx)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_k3_fd)
	temp1=vf(x,v1,v2,kshape)
	f1=deriv_copyldd(temp1,nx,dim=2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"gev_k3_f2fa=function(x,v1,v2,kshape){
	nx=length(x)
	f2=array(0,c(2,2,nx))

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_k3_fdd)
	temp1=vf(x,v1,v2,kshape)
	f2=deriv_copyldd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_k3_ldda=function(x,v1,v2,kshape){
	nx=length(x)
	ldd=matrix(0,2,2)

	kshape=movexiawayfromzero(kshape)

	vf=Vectorize(gev_k3_logfdd)
	temp1=vf(x,v1,v2,kshape)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_k3_lddda=function(x,v1,v2,kshape){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(gev_k3_logfddd)

	kshape=movexiawayfromzero(kshape)

	temp1=vf(x,v1,v2,kshape)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
