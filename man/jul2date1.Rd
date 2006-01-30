\name{jul2date1}
\alias{jul2date1}
\title{Converts Julian date into string date}
\description{
	Converts Julian day of year (DOY) into a string date 'DD.MM.YYYY'.
	If y is missing a non-leap year is assumed.
}
\usage{jul2date1(d,y)}
\arguments{
	\item{d}{DOY, numeric coerced into an integer.}
	\item{y}{Year, numeric coerced into an integer, default 2000.}
}
\value{
	\item{date}{Date, as character string 'DD.MM.YYYY'}	
}
\author{
	Joerg Schaber
}
\examples{
	jul2date1(151,1970)
}
\keyword{misc}
\keyword{utilities}
