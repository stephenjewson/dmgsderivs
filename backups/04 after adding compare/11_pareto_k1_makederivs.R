#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(extraDistr)
library(actuar)
#
# v1=shape (a)
# v2=scale (b)
f=function(x,v1,v2){v1*(v2^v1)/((x)^(v1+1))}
#f=function(x,v1,v2){v1*(v2^v1)/((x+v2)^(v1+1))} #this is the actuar version, which is different
compare("d",extraDistr::dpareto(4,3,2),f(4,3,2))
pareto_k1_fd=Deriv(f,"v1",nderiv=1)
pareto_k1_fdd=Deriv(f,"v1",nderiv=2)
#
logf=function(x,v1,v2){log(v1)+v1*log(v2)-(v1+1)*log(x)}
compare("l",extraDistr::dpareto(4,3,2,log=TRUE),logf(4,3,2))
pareto_k1_logfdd=Deriv(logf,"v1",nderiv=2)
pareto_k1_logfddd=Deriv(logf,"v1",nderiv=3)
#
sink("11c_pareto_k1_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k1_fd=")
print.function(pareto_k1_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("pareto_k1_fdd=")
print.function(pareto_k1_fdd)
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

