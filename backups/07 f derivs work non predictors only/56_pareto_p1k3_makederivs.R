#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
logf=function(x,t,v1,v2,v3){-v1-t*v2+exp(-v1-t*v2)*log(v3)-(exp(-v1-t*v2)+1)*log(x)}
pareto_p1k3_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
pareto_p1k3_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("56c_pareto_p1k3_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k3_logfdd=")
print.function(pareto_p1k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_p1k3_logfddd=")
print.function(pareto_p1k3_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_p1k3_ldda=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k3_logfdd)
	temp1=vf(x,t,v1,v2,kscale)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_p1k3_lddda=function(x,t,v1,v2,kscale){
	nx=length(x)
	vf=Vectorize(pareto_p1k3_logfddd)
	temp1=vf(x,t,v1,v2,kscale)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

