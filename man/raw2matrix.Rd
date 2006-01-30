\name{raw2matrix}
\alias{raw2matrix}
\title{Converts a numeric data frame to matrix}
\description{
	Converts a numeric data frame D with three columns (x, factor 1, factor 2)
	to a numeric matrix M where rows are ranks of levels of factor 1 and columns are
	ranks of levels of factor 2, missing values are set to NA.	
}
\usage{raw2matrix(D)}
\arguments{
	\item{D}{Data frame with three columns (x, factor 1, factor 2)}
}
\value{
	\item{M}{Numeric matrix where rows are ranks of levels of factor 1 and columns are
	ranks of levels of factor 2, missing values are set to NA.}
}
\author{
	Joerg Schaber
}
\examples{
	data(DWD)
	raw2matrix(DWD)
}
\keyword{utilities}
\keyword{misc}
