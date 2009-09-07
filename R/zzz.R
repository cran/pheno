.packageName <- "pheno"

.First.lib <- function(libname, pkgname) {
        library.dynam("pheno", pkgname, libname)
        print("pheno 1.5 library loaded")
	require(SparseM)
	require(quantreg)
	require(nlme)
}
