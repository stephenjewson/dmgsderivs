#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
# following my ordering here which is mu,sigma,df
logf=function(x,t,v1,v2,v3,v4){log(gamma(0.5*(v4+1)))-0.5*log(pi)-0.5*log(v4)-
		log(v3)-log(gamma(0.5*v4))-0.5*(v4+1)*log(1+(x-v1-t*v2)^2/(v4*v3*v3))}
lst_p1k4_logfdd=Deriv(logf,c("v1","v2","v3"),nderiv=2)
lst_p1k4_logfddd=Deriv(logf,c("v1","v2","v3"),nderiv=3)
#
sink("63c_lst_p1k4_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p1k4_logfdd=")
print.function(lst_p1k4_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_p1k4_logfddd=")
print.function(lst_p1k4_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p1k4_ldda=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	ldd=matrix(0,3,3)
	vf=Vectorize(lst_p1k4_logfdd)
	temp1=vf(x,t,v1,v2,v3,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_p1k4_lddda=function(x,t,v1,v2,v3,kdf){
	nx=length(x)
	lddd=array(0,c(3,3,3))
	vf=Vectorize(lst_p1k4_logfddd)
	temp1=vf(x,t,v1,v2,v3,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
