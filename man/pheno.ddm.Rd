\name{pheno.ddm}
\alias{pheno.ddm}
\title{Dense design matrix for phenological data}
\description{
	Creation of dense two-way classification design matrix
	for usage in robust parameter estimation with rq.fit.sfn (package nprq).
	The sum of the second factor is constrained to be zero. No general mean.
}
\usage{pheno.ddm(D,na.omit=TRUE)}
\arguments{
	\item{D}{Data frame with three columns: (observations, factor 1, factor 2).}
	\item{na.omit}{Determined whether missing values should be omitted or not. Default is TRUE.}
}
\details{
	In phenological applications observations should be the julian day
	of observation of a certain phase, factor 1 should be the observation year
	and factor 2 should be a station-id.
	Usually this is much easier created by:
	\code{y <- factor(f1)
	s <- factor(f2)
	ddm <- as.matrix.csr(model.matrix(~ y + s -1, contrasts=list(s=("contr.sum"))))}.
	However, this procedure can be quite memory demanding and might exceed storage
	capacity for large problems. 
	This procedure here is much less memory comsuming.
}
\value{
	\item{ddm}{Dense roworder matrix, matrix.csr format (see matrix.csr in package SparseM)}
	\item{D}{Data frame D sorted first by f2 then by f1 and with rows containing NA's removed.}
	\item{na.rows}{Rows in D that were omitted due to missing values.}
}
\author{
	Joerg Schaber
}
\seealso{
	\code{\link{model.matrix}}
	\code{\link[SparseM]{matrix.csr}}
}
\examples{
	data(DWD)
	ddm1 <- pheno.ddm(DWD)
	attach(DWD)
	y <- factor(DWD[[2]])
	s <- factor(DWD[[3]])
	ddm2 <- as.matrix.csr(model.matrix(~ y + s -1, contrasts=list(s=("contr.sum"))))
	identical(ddm1$ddm,ddm2)
}
\keyword{robust}
\keyword{design}
\keyword{models}
