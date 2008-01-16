\name{jul2date2}
\alias{jul2date2}
\title{Converts Julian date to integers day,month,year}
\description{
	Converts Julian day of year (DOY) into an integer date (day,month,year).
	If y is missing a non-leap year is assumed.
}
\usage{jul2date2(d,y)}
\arguments{
	\item{d}{DOY, numeric coerced into an integer.}
	\item{y}{Year, numeric coerced into an integer, default 2000.}
}
\value{
	\item{day}{Day of month as integer.}	
	\item{month}{Month of year as integer.}	
	\item{year}{Year as integer.}	
}
\author{
	Joerg Schaber
}
\examples{
	jul2date2(151,1970)
}
\keyword{misc}
\keyword{utilities}
