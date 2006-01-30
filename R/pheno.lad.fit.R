# Robust (least absolute deviations LAD/L1) fit of a two-way linear model 
# given a data frame with three columns (x, factor 1, factor 2) or 
# a matrix M where rows of M are ranks of factor 1 levels
# and columns of M are ranks of factor 2 levels, missing values are assumed to be NA or 0.
# No general mean and sum of f2 is constrained to be zero. 
# Estimation method: interior point method in case of dense implementation,
# else Barrodale-Roberts
# Depending on the size of the problem (n>1000) a dense matrix implementation is used.
# Some parameters of the dense algorithm are set according to my experience and should be OK
# for most cases. However, they are only tested for sparse design matrices up to the order
# ~ 90.000x2.900
# Output: p1 : parameter estimations of factor 1 (year effects)
# 		  p2 : parameter estimations of factor 2 (station effects)
#		  resid : residuals
#		  ierr : return code of the l1 estimation
#		  D  : the input as ordered data frame, ordered first after f2 then f1
pheno.lad.fit <- function(D) {
	if(!is.data.frame(D) && !is.matrix(D)) {
		stop("lad.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(D) && length(D)!=3) {
		stop("lad.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.matrix(D)) {
		D <- matrix2raw(D)
	}

	D <-  D[order(D[[3]],D[[2]]),]

	s <- factor(D[[3]])
	y <- factor(D[[2]])
	o <- as.vector(D[[1]],"numeric")
	nobs <- length(o)

	if(nobs > 1000) { # sparse matrix implementation
		S <- pheno.ddm(D)
		m <- S$ddm@dimension[2]
	    nnzdmax <- S$ddm@ia[nobs + 1] - 1
		l1fit <- rq.fit.sfn(S$ddm,o,tau=0.5,tmpmax=1000*m,nnzlmax=100*nnzdmax,small=1e-04)
		p1 <- as.vector(l1fit$coef[1:nlevels(y)],"numeric")
		p2 <- as.vector(contr.sum(nlevels(s)) %*% l1fit$coef[(nlevels(y)+1):(nlevels(y)+nlevels(s)-1)],"numeric")
		resid <- S$D[[1]] - l1fit$coef[match(y,levels(y))]
		ierr <- l1fit$ierr 
	}
	else {			# normal fit
		
		ddm <- as.matrix.csr(model.matrix(~ y + s -1, contrasts=list(s=("contr.sum"))))
		l1fit <- rq.fit(ddm,o,tau=0.5,method="br")
		p1 <-  as.vector(l1fit$coef[1:nlevels(y)],"numeric")
		p2 <- as.vector(contr.sum(nlevels(s)) %*% l1fit$coef[(nlevels(y)+1):(nlevels(y)+nlevels(s)-1)],"numeric")
		resid <- as.vector(residuals(l1fit),"numeric")
		ierr <- NA
	}
	return(list(p1=p1, p2=p2, resid=resid, ierr=ierr, D=D))
}
