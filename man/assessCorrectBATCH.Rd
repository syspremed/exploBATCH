\name{assessCorrectBATCH}
\alias{assessCorrectBATCH}
\title{
Assessing data corrected for batch effect using correctBATCH.
}
\description{
Assess if batch effect still exists after correctBATCH batch correction.
}
\usage{
assessCorrectBATCH(rerun2, qopt, nt, rerun22, theme)
}
\arguments{
  \item{rerun2}{
Results of running findBATCH to assess batch effect.
}
  \item{qopt}{
Optimal number of PCs.
}
  \item{nt}{
Number of biological classes.
}
  \item{rerun22}{
Results of running findBATCH to assess biological effect.
}
  \item{theme}{
ggplot theme
}
}
\details{
Assess if batch effect still exists after correctBATCH batch correction.
}
\value{
  Results in assessCorrectBATCH folder in the working directory.
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
