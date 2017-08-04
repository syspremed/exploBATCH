\name{ppcca}
\alias{ppcca}
\title{
Fitting PPCCA model
}
\description{
Fitting PPCCA model
}
\usage{
ppcca(Y, Yc, tYc, X, S, sumdiS, muhat, tp, tmp, minq = 2, maxq = 3, eps = 0.8, n.cores)
}
\arguments{
  \item{Y}{
Expression data
}
 \item{Yc}{
Centered and scaled data
}
  \item{tYc}{
Transposed centered and scaled data
}
  \item{X}{
Covariates
}
  \item{S}{
Sample covariance
}
  \item{sumdiS}{
Sum of the diagonal elements of sample covariance
}
  \item{muhat}{
Average expression of all features
}
  \item{tp}{
Number of dimensions from eigen
}
  \item{tmp}{
Eigen decomposion
}
  \item{minq}{
Minimum number of PCs considered
}
  \item{maxq}{
Maximum number of PCs considered
}
  \item{eps}{
Smallest value for convergence assessment
}
  \item{n.cores}{
Number of cores}
}
\details{
Fitting PPCCA model
}
\value{
Fitting PPCCA model
}
\references{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A., 2017. A novel and robust statistical method to diagnose and correct batch effects in genomic data.
}
\author{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line