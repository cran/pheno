\name{maxConnectedSet}
\alias{maxConnectedSet}
\title{Maximal connected set in a matrix}
\description{
	Finds connected data set, i.e. connected rows and columns
	of a numeric matrix M, that has the largest number of data entries.

}
\usage{maxConnectedSet(M)}
\arguments{
	\item{M}{Numeric matrix with missing values considered as 0, or
	a data frame. The data frame is internally converted to a matrix
   	and should have three columns (x, factor 1, factor 2) where
	x are considered the entries of the matrix, rows correspond to 
	levels of factor 2 and columns correspond to levels of factor 1.}
}
\details{
	In a two-way classification of linear models sometimes independent
	sets of normal equations are obtained due to missing data in the
	experiments design, i.e. the complete design matrix is not of full rank
	and thus no solution can be found. However, solutions of the independent
    sets of normal equations can still exist.
    This phenomenon is called 'connectedness' of the data.
    Especially in phenological analysis experimental designs are almost
    always unbalanced because of missing data. Thus, when combined time
    series are to be estimated, it is worth checking for and finding
    connected data sets for which combined time series can then be estimated.
	This can also be interpreted in the way that a prerequisite to obtain
	a combined time series is to have overlapping time series.
	Example (also see example data(Searle) from Searle (1997), page 324 and
	example in 'connectedSets'):
	In the following matrix dots represent missing values, X represent observations
	and the lines join the connected sets:

	:	X\_\_\_.\_\_\_.\_\_\_.\_\_\_X   .   .   . 		\cr
	:                       \eqn{\mid}{|}				\cr
	:	.   .   X\_\_\_.\_\_\_!\_\_\_.\_\_\_.\_\_\_X	\cr
	:                       \eqn{\mid}{|}           \eqn{\mid}{|}\cr
	:	.   X\_\_\_.\_\_\_.\_\_\_!\_\_\_X\_\_\_X   !\cr
	:                       \eqn{\mid}{|}       \eqn{\mid}{|}   \eqn{\mid}{|}\cr
	:	.   X\_\_\_.\_\_\_.\_\_\_!\_\_\_X\_\_\_X   !\cr
	:                       \eqn{\mid}{|}           \eqn{\mid}{|}\cr
	:	.   .   .   .   X   .   .   !					\cr
	:                       \eqn{\mid}{|}           \eqn{\mid}{|}\cr
	:	.   .   X\_\_\_.\_\_\_!\_\_\_.\_\_\_.\_\_\_X	\cr
	:                       \eqn{\mid}{|}				\cr
	:	.   .   .   X\_\_\_X   .   .   .				\cr

	Thus, in this matrix observations of rows 1, 5 and 7 or colums 1, 4 and 5 form
   	one connected set. Likewise observations of rows 2 and 6 (or columns 3 and 8)
   	and rows 3 and 4 (or columns 2, 6 and 7) form also connected sets, respectively.
}
\value{
	\item{ms}{maximal connected set as matrix or data frame, corresponding to the input.}
	\item{maxl}{Number of observations in the maximal connected data set.}
	\item{nsets}{Number of connected data sets.}
	\item{lsets}{Vector with number of observations in each connected data sets, i.e. lsets[i] is the number of observations in connected data set i.}
}
\references{
	Searle (1997) 'Linear Models'. Wiley. page 318.
}
\author{
	Joerg Schaber
}
\seealso{
	\code{\link{connectedSets}}
}
\examples{
	data(Searle)
	maxConnectedSet(Searle)
}
\keyword{design}
\keyword{models}
