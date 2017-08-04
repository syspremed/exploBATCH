\name{ppccaBeforeCorrection}
\alias{ppccaBeforeCorrection}
\title{
PPCCA before correction.
}
\description{
PPCCA before correction.
}
\usage{
ppccaBeforeCorrection(res1, grps, cGsub, batchCL, Conf, type, Ys, theme)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{res1}{
Results of PPCCA. 
}
  \item{grps}{
Number of batches.
}
\item{cGsub}{
  Colors
}
  \item{batchCL}{
A vector identifying batches.
}
  \item{Conf}{
A vector identifying biological classes.
}
  \item{type}{
Names of biological classes.
}
  \item{Ys}{
Expression data
}
  \item{theme}{
ggplot theme
}
}
\details{
PPCCA before correction to assess the optimal number of PCs using BIC.
}
\value{
  Results in ppccaBeforeCorrection folder in the working directory.
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