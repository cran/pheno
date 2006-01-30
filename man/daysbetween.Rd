\name{daysbetween}
\alias{daysbetween}
\title{Number of days between two dates}
\description{
	Number of days between date1 and date2.
}
\usage{daysbetween(d1,d2)}
\arguments{
	\item{d1}{Date as a character string 'DD.MM.YYYY'.}
	\item{d2}{Date as s character string 'DD.MM.YYYY'.}
}
\value{
	\item{ndays}{Number of days between d1 and d2.}	
}
\author{
	Joerg Schaber
}
\examples{
	daysbetween('31.05.1970','10.03.2004')
}
\keyword{misc}
\keyword{utilities}
