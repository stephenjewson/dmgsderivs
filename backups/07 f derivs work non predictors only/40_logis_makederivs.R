#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# note that the R definition of the pdf seems to be missing minus signs, so using wikipedia
#
f=function(x,v1,v2){(1/v2)*exp(-(x-v1)/v2)/((1+exp(-(x-v1)/v2))^2)}
compare("d",dlogis(1,2,3),f(1,2,3))
logis_fd=Deriv(f,c("v1","v2"),nderiv=1)
logis_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,v1,v2){1/(1+exp(-(x-v1)/v2))}
compare("p",plogis(1,2,3),p(1,2,3))
logis_pd=Deriv(p,c("v1","v2"),nderiv=1)
logis_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2){-log(v2)-(x-v1)/v2-2*log(1+exp(-(x-v1)/v2))}
compare("l",dlogis(1,2,3,log=TRUE),logf(1,2,3))
logis_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
logis_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("40c_logis_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("logis_fd=")
print.function(logis_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("logis_fdd=")
print.function(logis_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("logis_pd=")
print.function(logis_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("logis_pdd=")
print.function(logis_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("logis_logfdd=")
print.function(logis_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("logis_logfddd=")
print.function(logis_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"logis_f1fa=function(x,v1,v2){
	vf=Vectorize(logis_fd)
	f1=vf(x,v1,v2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"logis_f2fa=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyfdd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"logis_ldda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_logfdd)
	temp1=vf(x,v1,v2)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"logis_lddda=function(x,v1,v2){
	nx=length(x)
	vf=Vectorize(logis_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
