\name{tau}
\alias{tau}
\title{Kendall's normalized tau}
\description{
	Kendall's normalized tau for time series x
}
\usage{tau(x)}
\arguments{
	\item{x}{Numeric vector x.}
}
\details{
	Implicitly assumes a equidistant time series x. 
}
\value{
	\item{t}{Kendall's normalized tau.}
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
