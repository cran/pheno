\name{daylength}
\alias{daylength}
\title{Daylength at julian day i on latitude l}
\description{
	Calculates daylength [h] and declination angle delta [radians]
	on day i [julian day of year] for latitude l [degrees].	
}
\usage{daylength(i,l)}
\arguments{
	\item{i}{Integer as julian day of year (1-365)}	
	\item{l}{Float as latitude [degress]}
}
\value{
	\item{dl}{daylength [h]}
	\item{delta}{declination angle [degrees]}
}
\author{
	Joerg Schaber
}
\examples{
	daylength(as.integer(120),63)
}
\keyword{misc}
\keyword{utilities}
