\name{correctBatch}
\alias{correctBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch correction approach based on the PPCCA model.
}
\description{
Batch effect is corrected by subtracting its effect in the affected dimensions or PCs. 
Both corrected and uncorrected dimensions are then used to predict the data.
}
\usage{
correctBatch(res11, designX, Ys, batchCL, ngrps, grps, Conf, nt, cGsub, comres, theme, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res11}{
  Results of fitting the PPCCA model. 
}
  \item{designX}{
  Design matrix for the covariates
}
  \item{Ys}{
 Expression data
}
  \item{batchCL}{
  A vector identifying batches.
}
  \item{ngrps}{
  Number of batches
}
  \item{grps}{
Batches
}
  \item{Conf}{
  A vector identifying biological classes.
}
  \item{nt}{
  Number of biological classes
}
  \item{cGsub}{
  Colors
}
  \item{comres}{
  Results of ComBat correction
}
  \item{theme}{
ggplot theme
}
  \item{type}{
Plot
}
}
\details{
Batch effect is corrected by subtracting its effect only in the affected dimensions or PCs. 
Both the corrected and uncorrected dimensions are then used to predict the data.
}
\value{
  Results in correctBatch folder in the working directory.
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