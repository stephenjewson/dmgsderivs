#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# following my ordering here which is mu,sigma,df
f=function(x,t1,t2,v1,v2,v3,v4,v5){log(gamma(0.5*(v5+1)))-0.5*log(pi)-0.5*log(v5)-
		v3-t2*v4-log(gamma(0.5*v5))-0.5*(v5+1)*log(1+(x-v1-t1*v2)^2/(v5*exp(2*v3+2*t2*v4)))}
lst_p12k5_logfdd=Deriv(f,c("v1","v2","v3","v4"),nderiv=2)
lst_p12k5_logfddd=Deriv(f,c("v1","v2","v3","v4"),nderiv=3)
#
sink("81c_lst_p12k5_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p12k5_logfdd=")
print.function(lst_p12k5_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p12k5_logfddd=")
print.function(lst_p12k5_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p12k5_ldda=function(x,t1,t2,v1,v2,v3,v4,kdf){
	nx=length(x)
	ldd=matrix(0,4,4)
	vf=Vectorize(lst_p12k5_logfdd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p12k5_lddda=function(x,t1,t2,v1,v2,v3,v4,kdf){
	nx=length(x)
	lddd=array(0,c(4,4,4))
	vf=Vectorize(lst_p12k5_logfddd)
	temp1=vf(x,t1,t2,v1,v2,v3,v4,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
