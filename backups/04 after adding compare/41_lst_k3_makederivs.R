#
# make derivative codes for fitdistpu, one model at a time
#
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
library(extraDistr)
#
# following my ordering here which is mu,sigma,df
f=function(x,v1,v2,v3){(1/sqrt(v3*pi))*(1/gamma(v3/2))*gamma((v3+1)/2)*(1/v2)*(1+(x-v1)^2/(v3*v2*v2))^(-(v3+1)/2)}
compare("d",extraDistr::dlst(1,4,2,3),f(1,2,3,4))
lst_k3_fd=Deriv(f,c("v1","v2"),nderiv=1)
lst_k3_fdd=Deriv(f,c("v1","v2"),nderiv=2)
#
logf=function(x,v1,v2,v3){log(gamma(0.5*(v3+1)))-0.5*log(pi)-0.5*log(v3)-
		log(v2)-log(gamma(0.5*v3))-0.5*(v3+1)*log(1+(x-v1)^2/(v3*v2*v2))}
compare("l",extraDistr::dlst(1,4,2,3,log=TRUE),logf(1,2,3,4))
lst_k3_logfdd=Deriv(logf,c("v1","v2"),nderiv=2)
lst_k3_logfddd=Deriv(logf,c("v1","v2"),nderiv=3)
#
sink("41c_lst_k3_derivs.R")
#
cat("######################################################################\n")
cat("#' First derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_k3_fd=")
print.function(lst_k3_fd)
cat("######################################################################\n")
cat("#' Second derivative of the density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_k3_fdd=")
print.function(lst_k3_fdd)
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_k3_logfdd=")
print.function(lst_k3_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("lst_k3_logfddd=")
print.function(lst_k3_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_k3_ldda=function(x,v1,v2,kdf){
	nx=length(x)
	ldd=matrix(0,2,2)
	vf=Vectorize(lst_k3_logfdd)
	temp1=vf(x,v1,v2,kdf)
	ldd=deriv_copyldd(temp1,nx,dim=2)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"lst_k3_lddda=function(x,v1,v2,kdf){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	vf=Vectorize(lst_k3_logfddd)
	temp1=vf(x,v1,v2,kdf)
	lddd=deriv_copylddd(temp1,nx,dim=2)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
