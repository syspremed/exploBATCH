\name{assessComBat}
\alias{assessComBat}
\title{
Assessing ComBat corrected data
}
\description{
Assess if batch effect still exists after ComBat batch correction.
}
\usage{
assessComBat(rerun1, qopt, nt, rerun12, theme)
}
\arguments{
  \item{rerun1}{
Results of running findBATCH to assess batch effect.
}
  \item{qopt}{
Optimal number of PCs
}
  \item{nt}{
Number of biological classes.
}
  \item{rerun12}{
Results of running findBATCH to assess biological effect.
}
  \item{theme}{
ggplot theme
}
}
\details{
Assess if batch effect still exists after ComBat batch correction.
}
\value{
  Results in assessComBat folder in the working directory.
}
\references{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A., 2017. A novel and robust statistical method to diagnose and correct batch effects in genomic data.
}
\author{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
