\name{date2jul2}
\alias{date2jul2}
\title{Converts a date (day,month,year) to Julian date}
\description{
	Converts an integer date (day,month,year) into a Julian day of year (DOY).
	If y is missing, 2000 is assumed.
}
\usage{date2jul2(d,m,y)}
\arguments{
	\item{d}{Day of month, numeric coecerd into an integer.}
	\item{m}{Month of year, numeric coerced into an integer.}
	\item{y}{Year, numeric coerced into an integer, default 2000.}
}
\value{
	\item{doy}{Day of year as integer.}	
	\item{year}{Year as integer.}	
}
\author{
	Joerg Schaber
}
\examples{
	date2jul2(31,5,1970)
}
\keyword{misc}
\keyword{utilities}
