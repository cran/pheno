\name{pheno.lad.fit}
\alias{pheno.lad.fit}
\title{Fits a robust two-way linear model} 
\description{
	Fits a robust two-way linear model. 
	The model assumes both factors (f1 and f2) to be fixed.
	Errors are assumed to be i.i.d. No general mean and sum of
	f2 is constrained to be zero.
}
\usage{pheno.lad.fit(D,limit=1000)}
\arguments{
	\item{D}{Data frame with three columns (x, f1, f2) or a matrix 
	where rows are ranks of factor f1 levels and columns are ranks 
	of factor f2 levels and missing values are assumed to be NA or 0.}	
	\item{limit}{Integer that determines which algorithm to use (see Details).}
}
\details{
	The function minimizes the least absolute deviations (LAD or L1 norm)
	of the residuals of a two-way linear model.
	This function is basically a wrapper for the \code{rq.fit()} or \code{rq.fit.sfn()}
	functions of the \code{quantreg} package, respectively,
	adapted for the estimation of combined phenological time series. 
	Depending on the size of the problem \code{length(x)<=limit} 
	either the \code{rq.fit()} function using the Barrodale-Roberts algorithm is used or 
	(length(x)>1000) the corresponding dense matrix implementation with 
	\code{rq.fit.sfn()} using the Interior-Point method.
	In phenological applications, x should be the julian day
	of observation of a certain phase, factor f1 should be the observation year
	and factor f2 should be a station-id. 
	For efficiency reasons, the linear model is calcualted for treatment contrasts 
	and the constraint that the sum of f2 is zero, is adjusted afterwards.
	Note that the input data is sorted before fitting, such that subsequent
	analyses using the input data should be done using the sorted output data frame.
	
}
\value{
	\item{f1}{Estimated parameters of factor f1, in phenology this is precisely the combined time series.}
	\item{f1.lev}{Levels of f1. Should be the same order as f1.}
	\item{f2}{Estimated parameters of factor f2, in phenology these are precisely the station effects.}
	\item{f2.lev}{Levels of f2. Should be the same order as f2.}
	\item{resid}{Residuals}
	\item{ierr}{For length(x) > 1000 this is the return error code of \code{rq.fit.sfn()}}
	\item{D}{The input as ordered data frame, ordered first by f2 then by f1}
	\item{fit}{The fitted rq.fit model object.}
}
\author{
	Joerg Schaber
}
\references{
	Rousseeuw PJ, Leroy AM (1987) 'Robust estimation and outlier detection'. Wiley.
	Schaber J, Badeck F-W (2002) 'Evaluation of methods for the combination of phenological time series and outlier detection'. Tree Physiology 22:973-982
}
\seealso{
	\code{\link[quantreg]{rq.fit}}
	\code{\link[quantreg]{rq.fit.sfn}}
}
\examples{
	data(DWD)
	R <- pheno.lad.fit(DWD)				# robust parameter estimation
	plot(levels(factor(R$D[[2]])),R$p1,type="l")	# plot combined time series
	R$D[R$resid >= 30,]				# observation whose residuals
							# are > 30 days (outliers)
}
\keyword{design}
\keyword{models}
\keyword{ts}
