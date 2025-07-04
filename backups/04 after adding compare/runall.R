#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
library(Deriv)
library(fdrtool)
library(extraDistr)
library(actuar)
library(tictoc)
#
nonpredictors=TRUE
#
compare=function(t,a,b){
	cat(" ",t,":",a,"-",b,"=",round(a-b,8),"\n")
}
#
# 18
i=0
if(nonpredictors){
	i=i+1;cat(i,": exp:\n")
	source("./10_exp_makederivs.R")				;cat(" exp done\n\n")

	i=i+1;cat(i,": pareto:\n")
	source("./11_pareto_k1_makederivs.R")	;cat(" pareto done\n\n")

	i=i+1;cat(i,": halfnorm:\n")
	source("./20_halfnorm_makederivs.R")	;cat(" halfnorm done\n\n")

	i=i+1;cat(i,": norm:\n")
	source("./30_norm_makederivs.R")			;cat(" norm done\n\n")

	i=i+1;cat(i,": gnorm_k3:\n")
	source("./32_gnorm_k3_makederivs.R")	;cat(" gnorm done\n\n")

	i=i+1;cat(i,": lnorm:\n")
	source("./35_lnorm_makederivs.R")			;cat(" lnorm done\n\n")

	i=i+1;cat(i,": logis:\n")
	source("./40_logis_makederivs.R")			;cat(" logis done\n\n")

	i=i+1;cat(i,": lst_k3:\n")
	source("./41_lst_k3_makederivs.R")		;cat(" lst_k3 done\n\n")

	i=i+1;cat(i,": cauchy:\n")
	source("./42_cauchy_makederivs.R")		;cat(" cauchy done\n\n")

	i=i+1;cat(i,": gumbel:\n")
	source("./50_gumbel_makederivs.R")		;cat(" gumbel done\n\n")

	i=i+1;cat(i,": frechet:\n")
	source("./51_frechet_k1_makederivs.R");cat(" frechet_k1 done\n\n")

	i=i+1;cat(i,": weibull:\n")
	source("./52_weibull_makederivs.R")		;cat(" weibull done\n\n")

	i=i+1;cat(i,": gev_k3:\n")
	source("./53_gev_k3_makederivs.R")		;cat(" gev_k3 done\n\n")

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
predictors=FALSE
# 16
if(predictors){
	source("./55_exp_p1_makederivs.R")
	source("./56_pareto_p1k3_makederivs.R")
	source("./60_norm_p1_makederivs.R")
	source("./61_lnorm_p1_makederivs.R")
	source("./62_logis_p1_makederivs.R")
	source("./63_lst_p1k4_makederivs.R")
	source("./64_cauchy_p1_makederivs.R")
	source("./70_gumbel_p1_makederivs.R")
	source("./71_frechet_p2k1_makederivs.R")
	source("./73_weibull_p2_makederivs.R")
	source("./74_gev_p1k4_makederivs.R")
	source("./80_norm_p12_makederivs.R")
	source("./81_lst_p12k5_makederivs.R")

	cat(" gev_p1...(takes 15 secs)\n")
	tic()
	source("./150_gev_p1_makederivs.R")
	toc()

	cat(" gev_p12...(takes 47 secs and generates 471 lines of code)\n")
	tic()
	source("./151_gev_p12_makederivs.R")
	toc()
}

gev123=FALSE
if(gev123){
	cat(" gev_p123...(takes 488 secs and generates 562 lines of code with 479 separate variables)\n")
	tic()
	source("./152_gev_p123_makederivs.R")
	toc()
}

# total = 34
# noting that norm_dmgs and lnorm_dmgs are missing from this list
# because they use the same derivatives as norm and lnorm

rm(list=ls())
