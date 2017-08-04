\name{emQ}
\alias{emQ}
\title{
EM algorithm for PPCA model.
}
\description{
EM algorithm for PPCA model.
}
\usage{
emQ(q, Yc, tYc, S, tmp, tp, p, Ip, Np, Nvp, C2p, picon, IQ, penalty, eps, repV, V)
}
\arguments{
  \item{q}{
Number of PCs
}
  \item{Yc}{
Centered and scaled data
}
  \item{tYc}{
Transposed centered and scaled data
}
  \item{S}{
Sample covariance
}
  \item{tmp}{
Eigen decomposion
}
  \item{tp}{
Number of dimensions from eigen
}
  \item{p}{
Number of features
}
  \item{Ip}{
Identity matrix of p by p
}
  \item{Np}{
Some constant; n by p
}
  \item{Nvp}{
Some constant n by v by p
}
  \item{C2p}{
Some constant for the prior of the error variance
}
  \item{picon}{
Some constant for normal density computation
}
  \item{IQ}{
Identity matrices for different q values
}
  \item{penalty}{
penalty term for BIC
}
  \item{eps}{
Smallest value for convergence assessment
}
  \item{repV}{
Vector for maximum number of EM iteration in a vector
}
  \item{V}{
Maximum number of EM iterations.
}
}
\details{
EM algorithm for PPCA model.
}
\value{
EM algorithm for PPCA model.
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