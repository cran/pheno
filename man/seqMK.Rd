\name{seqMK}
\alias{seqMK}
\title{Sequential Mann-Kendall test for time series.}
\description{
	The sequential Mann-Kendall test on time series x detects 
	approximate potential trend turning points in time series.
}
\usage{seqMK(x)}
\arguments{
	\item{x}{Numeric vector x.}
}
\details{
	Implicitly assumes a equidistant time series x. 
	Calculates a progressive and a retrograde series of Kendall normalized tau's. 
	Points where the two lines cross are considered as approximate potential 
	trend turning points. When either the progressive or retrograde row exceed
	certain confidence limits before and after the crossing points, this trend
	turning point is considered significant at the corresponding level, 
	i.e. 1.96 for 95% significance and 2.58 for 99% significance etc.
}
\value{
	\item{prog}{Progressive row of Kendall's normalized tau's}
	\item{retr}{Retrograde row of Kendall's normalized tau's}
	\item{tp}{Boolean vector indicating at what indices of the original timeseries the prog and retr cross, i.e. TRUE at potential trend turning points.}
}
\author{
	Joerg Schaber
}
\references{
	Kendall M, Gibbons JD (1990) 'Rank correlation methods'. Arnold.
	Sneyers R (1990) 'On statistical analysis of series of observations. Technical  Note
		No 143. Geneva. Switzerland. World Meteorological Society.
	Schaber J (2003) 'Phenology in German in the 20th Century: Methods, analyses and models.
		Ph.D. Thesis. University of Potsdam. Germany.
		\url{http://pub.ub.uni-potsdam.de/2002meta/0022/door.htm}
}
\keyword{ts}
\keyword{utilities}
