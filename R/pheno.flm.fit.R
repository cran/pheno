# Fits a two-way linear model given a data frame with three columns
# (x, factor 1, factor 2) or a matrix M where rows of M are ranks of factor 1 levels
# and columns of M are ranks of factor 2 levels, missing values are assumed to be NA or 0.
# The model assumes the both factors f1 and f2 to be fixed.
# Errors are assumed to be i.i.d. No general mean and sum of f2 is constrained
# to be zero. For large problems, the sparse matrix implementation slm.fit is used.
# Output: f1 : coefficients of first factor
#        f1.lev: levels of first factor
#	 f2 :  coefficients of second fector
#        f2.lev: levels of second factor
#	 resid : residuals
#	 lclf1  : lower 95% confidence limits of factor f1
#	 uclf1  : upper 95% confidence limits of factor f1
#	 lclf2  : lower 95% confidence limits of factor f2
#	 uclf2  : upper 95% confidence limits of factor f2
#	 fit  : the fitted model object
pheno.flm.fit <- function(D) {
	if(!is.data.frame(D) && !is.matrix(D)) {
		stop("flm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(D) && length(D)!=3) {
		stop("flm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.matrix(D)) {
		D <- matrix2raw(D)
	}

	if(length(D[[1]]) > 1000) { # sparse matrix implementation
		S <- pheno.ddm(D)
		s <- factor(S$D[[3]])
		y <- factor(S$D[[2]])
		o <- as.vector(S$D[[1]],"numeric")
		nobs <- length(o)

		m <- S$ddm@dimension[2]
		fit <- slm.fit(S$ddm,o,tmpmax=1000*m,small=1e-06)

		fit$terms <- terms(o ~ y + s, contrasts=list(s=("contr.sum")),na.action=na.exclude)

		ny <- nlevels(y)
		f1 <- as.vector(summary.slm(fit)$coef[1:ny,1],"numeric")	
		f1.se <- as.vector(summary.slm(fit)$coef[1:ny,2],"numeric")	
		f2 <- as.vector(summary.slm(fit)$coef[-(1:ny),1],"numeric")	
		f2.se <- as.vector(summary.slm(fit)$coef[-(1:ny),2],"numeric")	
		resid <- as.vector(residuals(fit),"numeric")

		df <- summary.slm(fit)$df[2]
		lclf1 <- f1+qt(0.975,df)*summary.slm(fit)$coef[1:ny,2]
		uclf1 <- f1-qt(0.975,df)*summary.slm(fit)$coef[1:ny,2]

		lclf2 <- f2+qt(0.975,df)*summary.slm(fit)$coef[-(1:ny),2]
		uclf2 <- f2-qt(0.975,df)*summary.slm(fit)$coef[-(1:ny),2]
	}
	else {		# normal fit
		s <- factor(D[[3]])
		y <- factor(D[[2]])
		o <- as.vector(D[[1]],"numeric")
		nobs <- length(o)

		fit <- lm(o ~ y + s - 1, contrasts=list(s=("contr.sum")),na.action=na.exclude)

		yind <- grep("y",names(coef(fit)))
		sind <- grep("s",names(coef(fit)))
		
		f1 <- as.vector(summary(fit)$coef[yind,1],"numeric")	
		f1.se <- as.vector(summary(fit)$coef[yind,2],"numeric")	
		f2 <- as.vector(summary(fit)$coef[sind,1],"numeric")	
		f2.se <- as.vector(summary(fit)$coef[sind,2],"numeric")	
		resid <- as.vector(residuals(fit),"numeric")

		lclf1 <- as.vector(confint(fit)[yind,1],"numeric")
		uclf1 <- as.vector(confint(fit)[yind,2],"numeric")
		lclf2 <- as.vector(confint(fit)[sind,1],"numeric")
		uclf2 <- as.vector(confint(fit)[sind,2],"numeric")
	}

	return(list(f1=f1,f1.se=f1.se,f1.lev=levels(y),f2=f2,f2.se=f2.se,f2.lev=levels(s),resid=resid,lclf1=lclf1,uclf1=uclf1,lclf2=lclf2,uclf2=uclf2,fit=fit))
}
