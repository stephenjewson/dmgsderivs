#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(extraDistr)
#
# following my ordering here which is mu,sigma,df
f=function(x,t1,t2,v1,v2,v3,v4,v5){
	(1/sqrt(v5*pi))*(1/gamma(v5/2))*gamma((v5+1)/2)*(1/exp(v3+v4*t2))*(1+(x-v1-v2*t1)^2/(v5*exp(2*v3+2*v4*t2)))^(-(v5+1)/2)
}
compare("d",extraDistr::dlst(1,8,4+5*2,exp(.6+.7*.3)),f(1,2,.3,4,5,.6,.7,8))
lst_p12k3_fd=Deriv(f,c("v1","v2","v3","v4"),nderiv=1)
lst_p12k3_fdd=Deriv(f,c("v1","v2","v3","v4"),nderiv=2)
#
cat("  no cdf\n")
#
logf=function(x,t1,t2,v1,v2,v3,v4,v5){
	log(gamma(0.5*(v5+1)))-0.5*log(pi)-0.5*log(v5)-
		v3-v4*t2-log(gamma(0.5*v5))-0.5*(v5+1)*log(1+(x-v1-v2*t1)^2/(v5*exp(2*v3+2*v4*t2)))
}
compare("l",extraDistr::dlst(1,8,4+5*2,exp(.6+.7*.3),log=TRUE),logf(1,2,.3,4,5,.6,.7,8))
lst_p12k3_logfdd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=2)
lst_p12k3_logfddd=Deriv(logf,c("v1","v2","v3","v4"),nderiv=3)
#
sink("81c_lst_p12k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p12k3_fd=")
print.function(lst_p12k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p12k3_fdd=")
print.function(lst_p12k3_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p12k3_logfdd=")
print.function(lst_p12k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p12k3_logfddd=")
print.function(lst_p12k3_logfddd)
cat("############################################################\n")
#
cat("#' The first derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p12k3_f1fa=function(x,t,v1,v2,v3,v4,v5,kdf){
	vf=Vectorize(lst_p12k3_fd)
	f1=vf(x,t,v1,v2,v3,v4,v5,kdf)
	return(f1)
}\n"
)
cat("############################################################\n")
#
cat("#' The second derivative of the density\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p12k3_f2fa=function(x,t,v1,v2,v3,v4,v5,kdf){
	nx=length(x)
	vf=Vectorize(lst_p12k3_fdd)
	temp1=vf(x,t,v1,v2,v3,v4,v5,kdf)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}\n"
)
cat("############################################################\n")
#
###cat("#' The first derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"lst_p12k3_p1fa=function(x,t,v1,v2,v3,v4,v5,kdf){
###	vf=Vectorize(lst_p12k3_pd)
###	p1=vf(x,t,v1,v2,v3,v4,v5,kdf)
###	return(p1)
###}\n"
###)
cat("############################################################\n")
#
###cat("#' The second derivative of the cdf\n")
###cat("#' @inheritParams manf\n")
###cat(
###"lst_p12k3_p2fa=function(x,t,v1,v2,v3,v4,v5,kdf){
###	nx=length(x)
###	vf=Vectorize(lst_p12k3_pdd)
###	temp1=vf(x,t,v1,v2,v3,v4,v5,kdf)
###	p2=deriv_copyfdd(temp1,nx,dim=4)
###	return(p2)
###}\n"
###)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p12k3_ldda=function(x,t,v1,v2,v3,v4,v5,kdf){
	nx=length(x)
	vf=Vectorize(lst_p12k3_logfdd)
	temp1=vf(x,t,v1,v2,v3,v4,v5,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p12k3_lddda=function(x,t,v1,v2,v3,v4,v5,kdf){
	nx=length(x)
	vf=Vectorize(lst_p12k3_logfddd)
	temp1=vf(x,t,v1,v2,v3,v4,v5,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
