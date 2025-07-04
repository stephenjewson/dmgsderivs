#
# make derivative codes for fitdistpu, one model at a time
#
rm(list=ls())
setwd(paste(Sys.getenv('HOME'),'/03 pn/05 statistics/fitdistpu/makederivatives/',sep=""))
library(Deriv)
library(tictoc)
#
nonpredictors=FALSE

# 18
if(nonpredictors){
	source("./10_exp_makederivs.R")
	source("./11_pareto_k1_makederivs.R")
	source("./20_halfnorm_makederivs.R")
	source("./30_norm_makederivs.R")
	source("./32_gnorm_makederivs.R")
	source("./35_lnorm_makederivs.R")
	source("./40_logis_makederivs.R")
	source("./41_lst_k3_makederivs.R")
	source("./42_cauchy_makederivs.R")
	source("./50_gumbel_makederivs.R")
	source("./51_frechet_k1_makederivs.R")
	source("./52_weibull_makederivs.R")
	source("./53_gev_k3_makederivs.R")
	source("./100_gamma_makederivs.R")
	source("./101_invgamma_makederivs.R")
	source("./102_invgauss_makederivs.R")
	source("./110_gev_makederivs.R")
	source("./120_gpd_k1_makederivs.R")
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

gev123=TRUE
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
