#p=function(x){fdrtool::phalfnorm(x)}
#dp=Deriv(~p(x),"x",drule[["fdrtool::phalfnorm"]] <- alist(x=fdrtool::dhalfnorm(x)))
#print.function(dp)
#
# 1
#
cat("\ntest1 (works):\n")
sinpi=function(x){sin(pi*x)}
cospi=function(x){cos(pi*x)}
dsinpi=Deriv(~sinpi(x^2),"x",drule[["sinpi"]] <- alist(x=cospi(x)))
print.function(dsinpi)
#
# 2
#
cat("\ntest2:\n")
sinpi=function(x){sin(pi*x)}
cospi=function(x){cos(pi*x)}
drule[["sinpi"]] <- alist(x=cospi(x))
dsinpi=Deriv(~sinpi(x^2),"x",drule=drule)
print.function(dsinpi)
#
# 3
#
cat("\ntest3:\n")
phalfnorm=function(x){fdrtool::phalfnorm(x)}
dhalfnorm=function(x){fdrtool::dhalfnorm(x)}
drule[["phalfnorm"]] <- alist(x=dhalfnorm(x))
dp=Deriv(~phalfnorm(x^2),"x",drule=drule)
print.function(dp)
#
# 4
#
cat("\ntest4:\n")
phalfnorm=function(x){fdrtool::phalfnorm(x)}
dhalfnorm=function(x){fdrtool::dhalfnorm(x)}
drule[["phalfnorm"]] = alist(x=dhalfnorm(x))
dp=Deriv(~phalfnorm(a*x),c("x","a"),drule=drule)
print.function(dp)
#
# 5
#
cat("\ntest5\n")
p=function(x){fdrtool::phalfnorm(x)}
d=function(x){fdrtool::dhalfnorm(x)}
drule[["p"]] = alist(x=d(x))
dp=Deriv(~p(a*x),c("x","a"),drule=drule)
print.function(dp)
#
# 6
#
rm(list=ls())
cat("\ntest6\n")
p=function(x)x
d=function(x){(2/pi)*exp(-x*x/pi)}
drule[["p"]] = alist(x=d(x))
dp=Deriv(~p(x*v1),c("x","v1"),drule=drule)
print.function(dp)
#
# 7
#
rm(list=ls())
cat("\ntest7\n")
pmyhalfnorm=function(x)x
dmyhalfnorm=function(x){(2/pi)*exp(-x*x/pi)}
drule[["pmyhalfnorm"]] = alist(x=dmyhalfnorm(x))
halfnorm_pd=Deriv(~pmyhalfnorm(x*v1),c("v1"),nderiv=1,drule=drule)
print.function(halfnorm_pd)
#halfnorm_pdd=Deriv(~pmyhalfnorm(x*v1),c("v1"),nderiv=2,drule=drule)
#print.function(halfnorm_pdd)

