#
# make derivative codes for fitdistcp, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistcp/dmgsderivs/',sep=""))
library(Deriv)
library(fdrtool)
library(extraDistr)
library(actuar)
library(tictoc)
#
nonpredictors=TRUE
predictors=TRUE
gev123=TRUE
#
compare=function(t,a,b){
	cat(" ",t,":",a,"-",b,"=",round(a-b,8),"\n")
}
#
# 18
i=0

############## testing zone ###############
#source("./010_exp_makederivs.R")				;cat(" exp done\n\n")
#stop()
#source("./110_gev_makederivs.R")			;cat(" gev done\n\n")
#source("./011_pareto_k2_makederivs.R")	;cat(" pareto done\n\n")
#source("./020_halfnorm_makederivs.R")	;cat(" halfnorm done\n\n")
#source("./030_norm_makederivs.R")			;cat(" norm done\n\n")
#source("./030_norm_makederivs.R")			;cat(" norm done\n\n")
#	source("./035_lnorm_makederivs.R")			;cat(" lnorm done\n\n")
#source("./151_gev_p12_makederivs.R")
#source("./152_gev_p123_makederivs.R")
#rm(list=ls())
#cat("stopping deliberately\n")
#source("./042_cauchy_makederivs.R")		;cat(" cauchy done\n\n")
#stop()
###########################################


if(nonpredictors){
	i=i+1;cat(i,": exp:\n")
	source("./010_exp_makederivs.R")				;cat(" exp done\n\n")

	i=i+1;cat(i,": pareto:\n")
	source("./011_pareto_k2_makederivs.R")	;cat(" pareto done\n\n")

	i=i+1;cat(i,": halfnorm:\n")
	source("./020_halfnorm_makederivs.R")	;cat(" halfnorm done\n\n")

	i=i+1;cat(i,": norm:\n")
	source("./030_norm_makederivs.R")			;cat(" norm done\n\n")

	i=i+1;cat(i,": gnorm_k3:\n")
	source("./032_gnorm_k3_makederivs.R")	;cat(" gnorm done\n\n")

	i=i+1;cat(i,": lnorm:\n")
	source("./035_lnorm_makederivs.R")			;cat(" lnorm done\n\n")

	i=i+1;cat(i,": logis:\n")
	source("./040_logis_makederivs.R")			;cat(" logis done\n\n")

	i=i+1;cat(i,": lst_k3:\n")
	source("./041_lst_k3_makederivs.R")		;cat(" lst_k3 done\n\n")

	i=i+1;cat(i,": cauchy:\n")
	source("./042_cauchy_makederivs.R")		;cat(" cauchy done\n\n")

	i=i+1;cat(i,": gumbel:\n")
	source("./050_gumbel_makederivs.R")		;cat(" gumbel done\n\n")

	i=i+1;cat(i,": frechet:\n")
	source("./051_frechet_k1_makederivs.R");cat(" frechet_k1 done\n\n")

	i=i+1;cat(i,": weibull:\n")
	source("./052_weibull_makederivs.R")		;cat(" weibull done\n\n")

	i=i+1;cat(i,": gev_k3:\n")
	source("./053_gev_k3_makederivs.R")		;cat(" gev_k3 done\n\n")

	i=i+1;cat(i,": gamma:\n")
	source("./100_gamma_makederivs.R")		;cat(" gamma done\n\n")

	i=i+1;cat(i,": invgamm:\n")
	source("./101_invgamma_makederivs.R")	;cat(" invgamma done\n\n")

	i=i+1;cat(i,": invgauss:\n")
	source("./102_invgauss_makederivs.R")	;cat(" invgauss done\n\n")

	i=i+1;cat(i,": gev:\n")
	source("./110_gev_makederivs.R")			;cat(" gev done\n\n")

	i=i+1;cat(i,": gpd_k1:\n")
	source("./120_gpd_k1_makederivs.R")		;cat(" gpd_k1 done\n\n")

	i=i+1;cat(i,": gpd_k13:\n")
	source("./120_gpd_k13_makederivs.R")	;cat(" gpd_k13 done\n\n")
}
#
# 16
if(predictors){
	i=i+1;cat(i,": exp_p1:\n")
	source("./055_exp_p1_makederivs.R")

	i=i+1;cat(i,": pareto_p1k2:\n")
	source("./056_pareto_p1k2_makederivs.R")

	i=i+1;cat(i,": norm_p1:\n")
	source("./060_norm_p1_makederivs.R")

	i=i+1;cat(i,":lnorm_p1:\n")
	source("./061_lnorm_p1_makederivs.R")

	i=i+1;cat(i,": logis_p1:\n")
	source("./062_logis_p1_makederivs.R")

	i=i+1;cat(i,": lst_p1k3:\n")
	source("./063_lst_p1k3_makederivs.R")

	i=i+1;cat(i,": cauchy_p1:\n")
	source("./064_cauchy_p1_makederivs.R")

	i=i+1;cat(i,": gumbel_p1:\n")
	source("./070_gumbel_p1_makederivs.R")

	i=i+1;cat(i,": frechet_p2k1:\n")
	source("./071_frechet_p2k1_makederivs.R")

	i=i+1;cat(i,": weibull_p2:\n")
	source("./073_weibull_p2_makederivs.R")

	i=i+1;cat(i,": gev_p1k31:\n")
	source("./074_gev_p1k3_makederivs.R")

	i=i+1;cat(i,": gev_p12k3:\n")
	source("./074_gev_p12k3_makederivs.R")

#	source("./080_norm_p12_makederivs.R")
#	source("./081_lst_p12k3_makederivs.R")

	i=i+1;cat(i,": gev_p1:\n")
	cat(" gev_p1...(takes 15 secs)\n")
	tic()
	source("./150_gev_p1_makederivs.R")
	toc()

	i=i+1;cat(i,": gev_p12:\n")
	cat(" gev_p12...(takes 47 secs and generates 471 lines of code)\n")
	tic()
	source("./151_gev_p12_makederivs.R")
	toc()
}

if(gev123){
	i=i+1;cat(i,": gev_p123:\n")
	cat(" gev_p123...(takes 488 secs and generates 562 lines of code with 479 separate variables)\n")
	tic()
	source("./152_gev_p123_makederivs.R")
	toc()
}

# total = 34
# noting that norm_dmgs and lnorm_dmgs are missing from this list
# because they use the same derivatives as norm and lnorm

rm(list=ls())
