\name{connectedSets}
\alias{connectedSets}
\title{Connected sets in a matrix}
\description{
	Finds connected data sets, i.e. connected rows and columns
	of a numeric matrix M.
}
\usage{connectedSets(M)}
\arguments{
	\item{M}{Numeric matrix with missing values assumed to be NA or 0.}
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
	Example (also see example data(Simple) and example in 'maxConnectedSets'):
	In the following matrix dots represent missing values, X represent observations
	and the lines join the connected sets:
	   
	:	X\_\_\_X   .   .	\cr
	:	    \eqn{\mid}{|}	\cr
	:	X\_\_\_X   .   .	\cr
	:						\cr
	:	.   .   X\_\_\_X	\cr
	   
	Thus, in this matrix observations in rows 1 and 2 or colums 1 and 2 form one connected
   	set. Likewise row 3 (or columns 3 and 4) form also one connected set.	
}
\value{
	\item{rowclasses}{Vector of set numbers of rows of M (>=1). A value of \code{-1} indicates a row with missing data.}
	\item{colclasses}{Vector of set numbers of columns of M(>=1).  A value of \code{-1} indicates a column with missing data.}
}
\references{
	Searle (1997) 'Linear Models'. Wiley. page 318.
}
\author{
	Joerg Schaber
}
\seealso{
	\code{\link{maxConnectedSet}}
	\code{\link{getConnectedSets}}
}
\examples{
	data(Simple)
	connectedSets(Simple)
}
\keyword{design}
\keyword{models}
