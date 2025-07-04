#
# test 1
# -calls a function, for one model at a time, that prints the output to screen
# -tests d,p,q,r
#
rm(list=ls())
#setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives',sep=""))
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
f=function(x,v1){log(2)+log(v1)-log(pi)-x*x*v1*v1/pi}
fdd=Deriv(f,"v1",nderiv=2)
fddd=Deriv(f,"v1",nderiv=3)
#
sink("20c_halfnorm_derivs.R")
#
cat("#' Second derivative of log density\n")
cat("#' @inheritParams manf\n")
cat("halfnorm_fdd=")
print.function(fdd)
cat("#########################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"halfnorm_ldda=function(x,v1){
	nx=length(x)
	ldd=matrix(0,1,1)
	ldd[1,1]=sum(fdd(x,v1))/nx
	return(ldd)
}\n"
)
cat("#########################################################\n")
cat("#' Third derivative of log density\n")
cat("#' @inheritParams manf\n")
cat("halfnorm_fddd=")
print.function(fddd)
cat("#########################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"halfnorm_lddda=function(x,v1){
	nx=length(x)
	lddd=array(0,c(1,1,1))
	lddd[1,1,1]=sum(fddd(x,v1))/nx
	return(lddd)
}\n"
)
#
closeAllConnections()

