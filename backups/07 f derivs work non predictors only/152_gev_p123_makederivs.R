#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/97 MyRpackages/fitdistpu/R',sep=""))
library(Deriv)
#
logf=function(x,t1,t2,t3,v1,v2,v3,v4,v5,v6){-v3-t2*v4-(1+1/(v5+t3*v6))*log(1+(v5+t3*v6)*(x-v1-t1*v2)/exp(v3+t2*v4))-
		(1+(v5+t3*v6)*(x-v1-t*v2)/exp(v3+t2*v4))^(-1/(v5+t3*v6))}
gev_p123_logfdd=Deriv(logf,c("v1","v2","v3","v4","v5","v6"),nderiv=2)
gev_p123_logfddd=Deriv(logf,c("v1","v2","v3","v4","v5","v6"),nderiv=3)
#
sink("152c_gev_p123_derivs.R")
#
cat("############################################################\n")
cat("#' Second derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_logfdd=")
print.function(gev_p123_logfdd)
cat("############################################################\n")
cat("#' Third derivative of the log density\n")
cat("#' Created by Stephen Jewson\n")
cat("#' using Deriv() by Andrew Clausen and Serguei Sokol\n")
cat("#' @inheritParams manf\n")
cat("gev_p123_logfddd=")
print.function(gev_p123_logfddd)
cat("############################################################\n")
#
cat("#' The second derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_ldda=function(x,t,v1,v2,v3,v4,v5,v6){
	nx=length(x)
	ldd=matrix(0,6,6)

# are any of the xi values too near to zero?
# in this case, I can't adjust xi, only v5 and v6
# so I adjust v5, which will adjust xi
	minxi=10^-7
	for (i in 1:nx){
		xi=v5+t[i,3]*v6
		if(abs(xi)<minxi){
			if(xi>=0){
				v5=v5+minxi
			} else {
				v5=v5-minxi
			}
		}
	}

	vf=Vectorize(gev_p123_logfdd)
	t1=t[,1]
	t2=t[,2]
	t3=t[,3]
	temp1=vf(x,t1,t2,t3,v1,v2,v3,v4,v5,v6)
	ldd=deriv_copyldd(temp1,nx,dim=6)
	return(ldd)
}\n"
)
cat("############################################################\n")
#
cat("#' The third derivative of the normalized log-likelihood\n")
cat("#' @inheritParams manf\n")
cat(
"gev_p123_lddda=function(x,t,v1,v2,v3,v4,v5,v6){
	nx=length(x)
	lddd=array(0,c(6,6,6))

# are any of the xi values too near to zero?
# in this case, I can't adjust xi, only v5 and v6
# so I adjust v5, which will adjust xi
	minxi=10^-7
	for (i in 1:nx){
		xi=v5+t[i,3]*v6
		if(abs(xi)<minxi){
			if(xi>=0){
				v5=v5+minxi
			} else {
				v5=v5-minxi
			}
		}
	}

	vf=Vectorize(gev_p123_logfddd)
	t1=t[,1]
	t2=t[,2]
	t3=t[,3]
	temp1=vf(x,t1,t2,t3,v1,v2,v3,v4,v5,v6)
	lddd=deriv_copylddd(temp1,nx,dim=6)
	return(lddd)
}\n"
)
#
closeAllConnections()

setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
