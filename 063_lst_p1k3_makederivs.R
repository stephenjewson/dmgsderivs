#
# make derivative codes for fitdistcp, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistcp/R',sep=""))
library(Deriv)
library(extraDistr)
#
# following my ordering here which is mu,sigma,df
f=function(x,t,v1,v2,v3,v4){(1/sqrt(v4*pi))*(1/gamma(v4/2))*gamma((v4+1)/2)*(1/v3)*(1+(x-v1-v2*t)^2/(v4*v3*v3))^(-(v4+1)/2)}
compare("d",extraDistr::dlst(1,6,3+4*2,5),f(1,2,3,4,5,6))
lst_p1k3_fd=Deriv(f,c("v1","v2","v3"),nderiv=1)
lst_p1k3_fdd=Deriv(f,c("v1","v2","v3"),nderiv=2)
#
cat("  no cdf\n")
#
logf=function(x,t,v1,v2,v3,v4){log(gamma(0.5*(v4+1)))-0.5*log(pi)-0.5*log(v4)-
		log(v3)-log(gamma(0.5*v4))-0.5*(v4+1)*log(1+(x-v1-v2*t)^2/(v4*v3*v3))}
compare("l",extraDistr::dlst(1,6,3+4*2,5,log=TRUE),logf(1,2,3,4,5,6))
lst_p1k3_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
lst_p1k3_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("063c_lst_p1k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat("lst_p1k3_fd=")
print.function(lst_p1k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("lst_p1k3_fdd=")
print.function(lst_p1k3_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat("lst_p1k3_logfdd=")
print.function(lst_p1k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat("lst_p1k3_logfddd=")
print.function(lst_p1k3_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @returns Vector\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p1k3_f1fa=function(x,t,v1,v2,v3,kdf){
	vf=Vectorize(lst_p1k3_fd,\"x\")
	f1=vf(x,t,v1,v2,v3,kdf)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p1k3_f2fa=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	vf=Vectorize(lst_p1k3_fdd,\"x\")
	temp1=vf(x,t,v1,v2,v3,kdf)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}\n"
)
cat("############################################################\n")
#
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"lst_p1k3_p1fa=function(x,t,v1,v2,v3,kdf){
###	vf=Vectorize(lst_p1k3_pd,\"x\")
###	p1=vf(x,t,v1,v2,v3,kdf)
###	return(p1)
###}\n"
###)
cat("############################################################\n")
#
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"lst_p1k3_p2fa=function(x,t,v1,v2,v3,kdf){
###	nx=length(x)
###	vf=Vectorize(lst_p1k3_pdd,\"x\")
###	temp1=vf(x,t,v1,v2,v3,kdf)
###	p2=deriv_copyfdd(temp1,nx,dim=3)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @returns Matrix\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p1k3_ldda=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	vf=Vectorize(lst_p1k3_logfdd,\"x\")
	temp1=vf(x,t,v1,v2,v3,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @returns 3d array\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p1k3_lddda=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	vf=Vectorize(lst_p1k3_logfddd,\"x\")
	temp1=vf(x,t,v1,v2,v3,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
