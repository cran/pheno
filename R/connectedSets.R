.packageName <- "pheno"

.First.lib <- function(lib, pkg) {
        library.dynam("pheno", pkg, lib)
        print("pheno library loaded")
}
# Finds connected data sets of a numeric matrix M
# non-entries are considered 0
# Returns two vectors: 
# rowclasses[0..maxnr-1] : Class number of the respective rows
# colclasses[0..maxnc-1] : Class number of the respective cols
connectedSets <- function(M) {
	if(!is.matrix(M)) stop("connectedSets: first argument must be a matrix. Exiting ...")
	maxnr <- dim(M)[1]
	maxnc <- dim(M)[2]

	res <- .C("connectivity",M=as.vector(t(M),"numeric"),nrows=as.integer(maxnr),ncols=as.integer(maxnc),rowclasses=vector("integer",maxnr),colclasses=vector("integer",maxnc),PACKAGE="pheno")
	
	attach(res)
	return(rowclasses,colclasses)

}
