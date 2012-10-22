\name{date2jul1}
\alias{date2jul1}
\title{Converts string date to Julian date}
\description{
	Converts a string date "DD.MM.YYYY" into a Julian day
	of year (DOY).
}
\usage{date2jul1(d)}
\arguments{
	\item{d}{Date as charater string 'DD.MM.YYYY'.}
}
\value{
	\item{doy}{Day of year as integer.}	
	\item{year}{Year as integer.}	
}
\author{
	Joerg Schaber
}
\examples{
	date2jul1('31.05.1970')
}
\keyword{misc}
\keyword{utilities}
