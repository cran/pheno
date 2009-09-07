# calculates maximal daylength maxdl [h] at a certain latitude lat [degrees]
maxdaylength <- function(l) {

	res <- .C("maxdaylength",l=as.double(l),maxdl=double(1),PACKAGE="pheno")

	return(res$maxdl)
}
