# calculates daylength dl [h] and declination delta  [radians]
# on day i for latitude lat [degrees]
# declination: angle between sun rays and equatorial plane for 
# the whole earth (-23 degrees - + 23 degrees)
daylength <- function(i,l) {
		
	res <- .C("daylength",l=as.double(l),i=as.integer(i),dl=double(1),delta=double(1),PACKAGE="pheno")
	
	attach(res)
	return(dl,delta)

}
