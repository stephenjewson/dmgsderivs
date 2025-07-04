#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1,v2){log(v1)+v1*log(v2)-(v1+1)*log(v2)}
pareto_k1_logfdd=Deriv(f,"v1",nderiv=2)
pareto_k1_logfddd=Deriv(f,"v1",nderiv=3)
#
sink("11c_pareto_k1_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k1_logfdd=")
print.function(pareto_k1_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k1_logfddd=")
print.function(pareto_k1_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k1_ldda=function(x,v1,kscale){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(pareto_k1_logfdd)
	ldd[1,1]=sum(vf(x,v1,v2=kscale))/nx
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"pareto_k1_lddda=function(x,v1,kscale){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	temp=0
	vf=Vectorize(pareto_k1_logfddd)
	lddd[1,1,1]=sum(vf(x,v1,v2=kscale))/nx
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

