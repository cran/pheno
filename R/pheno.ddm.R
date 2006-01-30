# Automatic creation of dense two-way classification design matrix
# for usage of dense robust estimation with rq.fit.sfn (package qunatreg).
# The sum of the second factor is constrained to be zero. No general mean.
# Usually this is much easier created by:
# y <- factor(f1)
# s <- factor(f2)
# ddm <- as.matrix.csr(model.matrix(~ y + s -1, contrasts=list(s=("contr.sum"))))
# however, this procedure is quite memory demanding and might exceed storage
# capacity for large problems. 
# This procedure here is much less memory comsuming.
# input: data frame with three columns: (observations, factor 1, factor 2)
# (phenology: observation day of phase, year, station)
# output: dense roworder matrix, matrix.csr format (see matrix.csr in package SparseM)
# and the sorted data frame D (data frame is being sorted first by f2 then by f1 )
pheno.ddm <- function(D) {
	if(!is.data.frame(D) || length(D) != 3) 
		stop("pheno.ddm: argument must be data frame with 3 fields. Exiting ...")

	# order first by factor 2 then by factor 1
	D <-  D[order(D[[3]],D[[2]]),]
	no <- length(D[[1]])	 	# number of observations
	f1 <- factor(D[[2]]) 		# factor 1: year
	n1 <- nlevels(f1) 			# number of levels factor 1 (phenology: years)
	f2 <- factor(D[[3]])		# factor 2: station
	n2 <- nlevels(f2)			# number of levels factor 2 (phenology: station)
	
	# ra: Object of class numeric, a real array of nnz elements containing the 
	#	non-zero elements of A, stored in row order. Thus, if i<j, all elements 
	#	of row i precede elements from row j. The order of elements within the 
	#	rows is immaterial. 
	# ja: Object of class integer, an integer array of nnz elements containing 
	#	the column indices of the elements stored in  ra. 
	# ia: Object of class integer, an integer array of n+1 elements containing 
	#	pointers to the beginning of each row in the arrays  ra  and  ja . 
	#	Thus  ia[i]  indicates the position in the arrays  ra  and  ja  where 
	#	the ith row begins. The last, (n+1)st, element of  ia  indicates 
	# 	where the n+1 row would start, if it existed. 
	# dimension: Object of class integer, dimension of the matrix
	
	# number of observations in last level of factor 2
	nols <- length(split(D,D[[3]])[[n2]][[1]]) 
	# number of non-zero elements in ra and ja 
	nnz <- (no-nols)*2+nols*n2
	ra=numeric(nnz)
	ja=integer(nnz)
	ia=integer(no+1)
	dimension=integer(2)
	ra[1:((no-nols)*2)] <- 1
	ra[((no-nols)*2+1):nnz] <- rep(c(1,rep(-1,n2-1)),nols)

	ia[1] <- as.integer(1) 
	for(i in 1:(no-nols)) {
		ja[2*i-1] <- as.integer(f1[[i]])
		ja[2*i] <- n1 + as.integer(f2[[i]])
		ia[i+1] <- as.integer(ia[i] + 2)
	}
	k <- 2*(no-nols)+1
	for(i in (no-nols+1):no) {
		ia[i+1] <- as.integer(ia[i] + n2)
		ja[k] <- as.integer(f1[[i]])
		ja[(k+1):(k+n2-1)] <- as.integer(c((n1+1):(n1+n2-1)))
		k <- k + n2
	}
	dim <- as.integer(c(no,n1+n2-1))
	ddm <- new("matrix.csr",ra=ra,ja=ja,ia=ia,dimension=dim)
	return(list(ddm=ddm,D=D))
}
