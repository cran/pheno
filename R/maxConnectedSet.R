# Returns largest connected data set of a either numeric data frame D 
# with three columns (x, factor 1, factor 2) or a n*m matrix M,
# where the n rows correspond to n levels of factor 2 and m columns
# correspond to m levels of factor 1.
# Output as data frame or matrix, depending on input, with number
# of data entries in the maximal connected set, and number of connectes sets
maxConnectedSet <- function(M) {
	if(!is.data.frame(M) && !is.matrix(M)) {
		stop("maxConnectedSet: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(M) && length(M)!=3) {
		stop("maxConnectedSet: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(M)) {
		f1 <- factor(M[[2]])
		f2 <- factor(M[[3]])
		M <- raw2matrix(M)
		out <- 1
	}
	else { out <- 0 }

	sets <- connectedSets(M)	# find connected sets

	nsets <- length(unique(sets$colclasses)) # number of connected sets

	lsets <- vector("numeric",nsets)

	# find sets with maximal numbers of data entries
	maxl <- 0
	for(i in unique(sets$colclasses)) {
		lsets[i] <- length(which(M[which(sets$rowclasses==i),which(sets$colclasses==i)]!=0))
		if(lsets[i] > maxl) { max <- i; maxl <- lsets[i] }
	}
	
	ms <- M[which(sets$rowclasses==max),which(sets$colclasses==max)]
	
	if(out == 0) { # return matrix
		return(ms,maxl,nsets,lsets)
	}
	else {			# return data frame
		ms <- matrix2raw(ms,as.vector(levels(f1),"numeric"), as.vector(levels(f2),"numeric"))
		return(ms,maxl,nsets,lsets)
	}
}
