############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_logfdd=function (x, v1, v2) 
{
    .e1 <- v2^2
    .e2 <- (2 * .e1)^2
    .e3 <- x - v1
    .e4 <- 1/.e1
    c(v1 = c(v1 = -.e4, v2 = -(8 * (v2 * .e3/.e2))), v2 = c(v1 = -(2 * 
        (.e3/v2^3)), v2 = .e4 + 4 * ((1 - 16 * (v2^4/.e2)) * 
        .e3^2/.e2)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @inheritParams manf
norm_logfddd=function (x, v1, v2) 
{
    .e1 <- 2 * v2^2
    .e2 <- .e1^2
    .e3 <- v2^4
    .e4 <- 16 * (.e3/.e2)
    .e5 <- v2^3
    .e6 <- x - v1
    .e7 <- 1 - .e4
    .e8 <- 2/.e5
    .e10 <- c(v1 = .e8, v2 = -(8 * (.e7 * .e6/.e2)))
    c(v1 = c(v1 = c(v1 = 0, v2 = 8 * (v2/.e2)), v2 = .e10), v2 = c(v1 = .e10, 
        v2 = c(v1 = 6 * (.e6/.e3), v2 = -(.e8 + 4 * (.e5 * (16 * 
            .e7 + 16 * (4 - .e4)) * .e6^2/.e1^4)))))
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @inheritParams manf
norm_ldda=function(x,v1,v2){
	nx=length(x)
	ldd=matrix(0,2,2)
	temp=0
	for (i in 1:nx){
		temp=temp+norm_logfdd(x[i],v1,v2)
	}
	for (i in 1:2){
		for (j in 1:2){
			ldd[i,j]=temp[(i-1)*2+j]/nx
		}
	}
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @inheritParams manf
norm_lddda=function(x,v1,v2){
	nx=length(x)
	lddd=array(0,c(2,2,2))
	temp=0
	for (i in 1:nx){
		temp=temp+norm_logfddd(x[i],v1,v2)
	}
	for (i in 1:2){
		for (j in 1:2){
			for (k in 1:2){
				lddd[i,j,k]=temp[(i-1)*4+(j-1)*2+k]/nx
			}
		}
	}
	return(lddd)
}
