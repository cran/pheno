\name{pheno.mlm.fit}
\alias{pheno.mlm.fit}
\title{Fits a two-way linear mixed model} 
\description{
	Fits a two-way linear mixed model. 
	The model assumes the first factor f1 to be fixed and the second factor f2 to
	be random. Errors are assumed to be i.i.d. No general mean and sum of 
	f2 is constrained to be zero.
}
\usage{pheno.mlm.fit(D)}
\arguments{
	\item{D}{Data frame with three columns (x, f1, f2) or a matrix
	where rows are ranks of factor f1 levels and columns are ranks 
	of factor f2 levels and missing values are set to 0.}	
}
\details{
	This function is basically a wrapper for the \code{lme()} function of
	the \code{nlme} package, adapted for the estimation of combined
	phenological time series. Estimation method: restricted maximum likelihood (REML)
	In phenological application, x should be the julian day
	of observation of a certain phase, factor f1 should be the observation year
	and factor f2 should be a station-id. 
}
\value{
	\item{fixed}{Estimated fixed effects, in phenology this is precisely the combined time series.}
	\item{fixed.lev}{Levels of fixed effects. Should be the same order as fixed effects.}
	\item{random}{Estimated random effects, in phenology these are the station effects.}
	\item{random.lev}{Levels of random effects. Should be the same order as random effects.}
	\item{SEf1}{Standard error group f1, i.e. square root of variance component fixed effect.}
	\item{SEf2}{Standard error group f2, i.e. square root of variance component random effect.}
	\item{lclf}{Lower 95 percent confidence limit of fixed effects.}
	\item{uclf}{Upper 95 percent confidence limit of fixed effects.}
	\item{fit}{The fitted lme model object.}
}
\references{
	Searle (1997) 'Linear Models'. Wiley.
	Schaber J, Badeck F-W (2002) 'Evaluation of methods for the combination of phenological time series and outlier detection'. Tree Physiology 22:973-982
}
\author{
	Joerg Schaber
}
\seealso{
	\code{\link[nlme]{lme}}
}
\examples{
	data(DWD)
	R <- pheno.mlm.fit(DWD)								# parameter estimation
	plot(levels(factor(DWD[[2]])),R$fixed,type="l")		# plot combined time series
	tr <- lm(R$fixed~rank(levels(factor(DWD[[2]]))))	# trend estimation
	summary(tr)$coef[2]									# slope of trend
	summary(tr)$coef[4]									# standard error of trend
}
\keyword{design}
\keyword{models}
\keyword{ts}
