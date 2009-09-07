\name{pheno.flm.fit}
\alias{pheno.flm.fit}
\title{Fits a two-way linear fixed model} 
\description{
	Fits a two-way linear fixed model. 
	The model assumes the first factor f1 the second factor f2 to be fixed. 
	Errors are assumed to be i.i.d. No general mean and sum of 
	f2 is constrained to be zero.
}
\usage{pheno.flm.fit(D,limit=1000)}
\arguments{
	\item{D}{Data frame with three columns (x, f1, f2) or a matrix
	where rows are ranks of factor f1 levels and columns are ranks 
	of factor f2 levels and missing values are assumed to be NA or 0.}	
	\item{limit}{Integer that determines which algorithm to use (see Details).}
}
\details{
	This function is basically a wrapper for the \code{slm.fit()} function
	form the \code{SparseM} package, adapted for the estimation of combined phenological time series.
	In phenological application, x should be the julian day
	of observation of a certain phase, factor f1 should be the observation year
	and factor f2 should be a station-id.
	For large problems \code{length(x)>limit}, the linear model is calculated
	for treatment contrasts for efficiency reasons, and the constraint that the sum of f2 is zero,
	is adjusted afterwards. This results in a slight over-estimation of
	standard errors.
	Note that the input data is sorted before fitting, such that subsequent
	analyses using the input data should be done using the sorted output data frame.
}
\value{
	\item{f1}{Estimated fixed effects f1, in phenology this is precisely the combined time series.}
	\item{f1.se}{f1 estimated standard error.}
	\item{f1.lev}{Levels of f1. Should be the same order as f1.}
	\item{f2}{Estimated fixed effects f2, in phenology these are the station effects.}
	\item{f2.se}{f2 estimated standard error.}
	\item{f2.lev}{Levels of f2. Should be the same order as f2.}
	\item{resid}{Residuals}
	\item{lclf1}{Lower 95 percent confidence limit of factor f1.}
	\item{uclf1}{Upper 95 percent confidence limit of factor f1.}
	\item{lclf2}{Lower 95 percent confidence limit of factor f2.}
	\item{uclf2}{Upper 95 percent confidence limit of factor f2.}
	\item{D}{The input as ordered data frame, ordered first by f2 then by f1}
	\item{fit}{The fitted lm model object.}
}
\references{
	Searle (1997) 'Linear Models'. Wiley.
	Schaber J, Badeck F-W (2002) 'Evaluation of methods for the combination of phenological time series and outlier detection'. Tree Physiology 22:973-982
}
\author{
	Joerg Schaber
}
\seealso{
	\code{\link[stats]{lm}}
}
\examples{
	data(DWD)
	R <- pheno.flm.fit(DWD)					# parameter estimation
}
\keyword{design}
\keyword{models}
\keyword{ts}
