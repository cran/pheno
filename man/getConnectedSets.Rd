\name{getConnectedSets}
\alias{getConnectedSets}
\title{Finds connected sets in a matrix or data frame}
\description{
	Finds a list of connected data sets in a matrix or data frame
	and returns them accordingly.

}
\usage{getConnectedSets(M)}
\arguments{
	\item{M}{Numeric matrix with missing values considered as 0, or
	a data frame. The data frame is internally converted to a matrix
   	and should have three columns (x, factor 1, factor 2) where
	x are considered the entries of the matrix, rows correspond to 
	levels of factor 2 and columns correspond to levels of factor 1.}
}
\details{
	getConnctedSets returns a list of connected data sets as numeric data frames D
	with three columns (x, factor 1, factor 2) or a n*m matrix M,
	where the n rows correspond to n levels of factor 2 and m columns
	correspond to m levels of factor the respective factors.
	Output as data frame or matrix, depending on input.
}
\value{
	\item{cs_i}{List of connected sets as matrix or data frame, corresponding to the input.
		     named as cs_i with i being the number of the connected sets.}
}
\references{
	Searle (1997) 'Linear Models'. Wiley. page 318.
}
\author{
	Joerg Schaber
}
\seealso{
	\code{\link{connectedSets}}
	\code{\link{maxConnectedSet}}
}
\examples{
	data(Searle)
	getConnectedSets(Searle)
}
\keyword{design}
\keyword{models}
