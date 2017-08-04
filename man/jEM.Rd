\name{jEM}
\alias{jEM}
\title{
EM algorithm for PPCCA model jackkniifing.
}
\description{
EM algorithm for PPCCA model jackkniifing.
}
\usage{
jEM(n, Y, idx, nj, p, Covars, Alphaopt, Wopt, Sigopt, Iq, repV, eps, Np, Nvp, C2p, Ip, picon, V)
}
\arguments{
  \item{n}{
Number of samples
}
  \item{Y}{
Expression data
}
  \item{idx}{
Identify samples
}
  \item{nj}{
Number of jackknifing samples
}
  \item{p}{
Number of features
}
  \item{Covars}{
Covariates
}
  \item{Alphaopt}{
Regression coefficients
}
  \item{Wopt}{
Loadings
}
  \item{Sigopt}{
Error variance
}
  \item{Iq}{
Identity matrix q by q}
  \item{repV}{
Vector for maximum number of EM iteration in a vector
}
  \item{eps}{
Smallest value for convergence assessment
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
  \item{Ip}{
Identity matrix of p by p
}
  \item{picon}{
Some constant for normal density computation
}
  \item{V}{
Maximum number of EM iterations.
}
}
\details{
EM algorithm for PPCCA model jackkniifing.
}
\value{
EM algorithm for PPCCA model jackkniifing.
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