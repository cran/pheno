\name{matrix2raw}
\alias{matrix2raw}
\title{Converts numeric matrix to data frame}
\description{
	Converts a numeric matrix M into a dataframe D
	with three columns (x, factor 1, factor 2)
	where rows of M are ranks of factor 1 levels and columns of M are
	ranks of factor 2 levels, missing values are set to 0.
}
\usage{matrix2raw(M,l1,l2)}
\arguments{
	\item{M}{Numeric matrix}
	\item{l1}{\bold{Optional} numeric vector of level names of column 2 (factor 1)
		of returned data frame. If missing it is assigned row numbers of \code{M}.}
	\item{l2}{\bold{Optional} numeric vector of level names of column 3 (factor 2)
		of returned data frame. If missing it is assigned column numbers of \code{M}.}
}
\value{
	\item{D}{Data frame with three columns: (y,f1,f1). \code{y}: observations, i.e.
	non-zero entries, in matrix. \code{f1}: factor 1, i.e. row number of M or \code{l1}.
	\code{f2}: factor 2, i.e. column number of M or \code{l2}. D is ordered first by
	factor 2 and then factor 1.}
}
\author{
	Joerg Schaber
}
\examples{
	data(DWD)
	M <- raw2matrix(DWD)	# conversion to matrix
	D1 <- matrix2raw(M)	# back conversion, but with different level names
	D2 <- matrix2raw(M,c(1951:1998),c(1:9))	# with original level names
}
\keyword{utilities}
\keyword{misc}
