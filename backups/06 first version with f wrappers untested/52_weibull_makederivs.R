#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# v1=shape
# v2=scale
f=function(x,v1,v2){(v1/v2)*((x/v2)^(v1-1))*exp(-((x/v2)^v1))}
compare("d",dweibull(1,2,3),f(1,2,3))
weibull_fd=Deriv(f,c("v1","v2"),nderiv=1)
weibull_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
p=function(x,v1,v2){1-exp(-((x/v2)^v1))}
compare("p",pweibull(1,2,3),p(1,2,3))
weibull_pd=Deriv(p,c("v1","v2"),nderiv=1)
weibull_pdd=Deriv(p,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2){log(v1)-log(v2)+(v1-1)*log(x/v2)-(x/v2)^v1}
compare("l",dweibull(1,2,3,log=TRUE),logf(1,2,3))
weibull_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
weibull_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("52c_weibull_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_fd=")
print.function(weibull_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_fdd=")
print.function(weibull_fdd)
cat("######################################################################\n")
cat("#' First derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_pd=")
print.function(weibull_pd)
cat("######################################################################\n")
cat("#' Second derivative of the cdf\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_pdd=")
print.function(weibull_pdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_logfdd=")
print.function(weibull_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("weibull_logfddd=")
print.function(weibull_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_f1fa=function(x,v1,v2){
	nx=length(x)
	f1=matrix(0,2,nx)
	vf=Vectorize(weibull_fd)
	temp1=vf(x,v1,v2)
	f1=deriv_copyldd(temp1,nx,dim=2)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_f2fa=function(x,v1,v2){
	nx=length(x)
	f2=array(0,c(2,2,nx))
	vf=Vectorize(weibull_fdd)
	temp1=vf(x,v1,v2)
	f2=deriv_copyldd(temp1,nx,dim=2)
	return(f2)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"weibull_ldda=function(x,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(weibull_logfdd)
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
"weibull_lddda=function(x,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(weibull_logfddd)
	temp1=vf(x,v1,v2)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
