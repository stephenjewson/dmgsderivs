#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1){v1*exp(-v1*x)}
compare("d",dexp(1,2),f(1,2))
exp_fd=Deriv(f,"v1",nderiv=1)
exp_fdd=Deriv(f,"v1",nderiv=2)
#
p=function(x,v1){1-exp(-v1*x)}
compare("p",pexp(1,2),p(1,2))
exp_pd=Deriv(p,"v1",nderiv=1)
exp_pdd=Deriv(p,"v1",nderiv=2)
#
logf=function(x,v1){log(v1)-v1*x}
compare("l",dexp(1,2,log=TRUE),logf(1,2))
exp_logfdd=Deriv(logf,"v1",nderiv=2)
exp_logfddd=Deriv(logf,"v1",nderiv=3)
#
sink("10c_exp_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_fd=")
print.function(exp_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_fdd=")
print.function(exp_fdd)
cat("######################################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_logfdd=")
print.function(exp_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("exp_logfddd=")
print.function(exp_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"exp_ldda=function(x,v1){
	nx=length(x)
	ldd=matrix(0,1,1)
	vf=Vectorize(exp_logfdd)
	ldd[1,1]=sum(vf(x,v1))/nx
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"exp_lddda=function(x,v1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	vf=Vectorize(exp_logfddd)
	lddd[1,1,1]=sum(vf(x,v1))/nx
	return(lddd)
}\n"
)
#
closeAllConnections()
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))

