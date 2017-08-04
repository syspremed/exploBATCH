\name{correctComBat}
\alias{correctComBat}
\title{
Batch correction with ComBat
}
\description{
Batch correction with ComBat
}
\usage{
correctComBat(Ys, batchCL, grps, cGsub, Conf, theme, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ys}{
  Expression data
}
  \item{batchCL}{
  A vector identifying batches.
}
  \item{grps}{
  Number of batches
}
  \item{cGsub}{
  Colors
}
  \item{Conf}{
  A vector identifying biological classes.
}
  \item{theme}{
ggplot theme
}
  \item{type}{
plot
}
}
\details{
Batch correction with ComBat
}
\value{
  Results in correctComBat folder in the working directory.
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