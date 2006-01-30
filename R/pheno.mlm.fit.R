# Fits a two-way linear mixed model given a data frame with three columns
# (x, factor 1, factor 2) or a matrix M where rows of M are ranks of factor 1 levels
# and columns of M are ranks of factor 2 levels, missing values are assumed to be NA or 0.
# The model assumes the first factor f1 to be fixed and the second factor f2 to
# be random. Errors are assumed to be i.i.d. No general mean and sum of f2 is constrained
# to be zero. Estimation method: restricted maximum likelihood (REML)
# Output: fixed : fixed effects
#		  random: random effects
#		  resid : residuals
#		  SEf1  : standard error group f1 (square root of variance component fixed effect)
#		  SEf2  : standard error group f2 (square root of variance component random effect)
#		  lclf1  : lower 95% confidence limit of fixed effects (factor f1)
#		  uclf1  : upper 95% confidence limit of fixed effects (factor f1)
#		  lclf2  : lower 95% confidence limit of random effects (factor f2)
#		  uclf2  : upper 95% confidence limit of ramdom effects (factor f2)
#		  n      : number of used data points
#		  df     : degrees of freedom
pheno.mlm.fit <- function(D) {
	if(!is.data.frame(D) && !is.matrix(D)) {
		stop("mlm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(D) && length(D)!=3) {
		stop("mlm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.matrix(D)) {
		D <- matrix2raw(D)
	}

	s <- factor(D[[3]])
	y <- factor(D[[2]])
	o <- as.vector(D[[1]],"numeric")
	remlfit <- lme(o ~ y - 1 ,random = ~ 1 | s, method="REML", contrasts=list(s=("contr.sum")),na.action=na.exclude)

	fixed <- as.vector(fixed.effects(remlfit),"numeric")
	random <- as.vector(random.effects(remlfit)[[1]],"numeric")
	resid <- as.vector(residuals(remlfit),"numeric")
	SEf1 <-  summary(remlfit)$sigma
	SEf2 <-  attr(remlfit$apVar,"Pars")[[2]]
	lclf1 <- as.vector(intervals(remlfit)$fixed[,1],"numeric")
	uclf1 <- as.vector(intervals(remlfit)$fixed[,3],"numeric")
	return(list(fixed=fixed, random=random, resid=resid, SEf1=SEf1, SEf2=SEf2, lclf1=lclf1, uclf1=uclf1, n=remlfit$dims$N, df=as.numeric(remlfit$fixDF$terms)))
}
