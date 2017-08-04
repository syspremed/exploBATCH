\name{findBATCH}
\alias{findBATCH}
\title{
A function to formally tests presence of batch effect.
}
\description{
It searches for directions in the principal subspace that are associated with the batch variable.
}
\usage{
findBATCH(res2, qopt, nt, rerun, theme)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res2}{
Results of running findBATCH to assess batch effect.
}
  \item{qopt}{
Optimal number of PCs
}
  \item{nt}{
Number of biological classes.
}
  \item{rerun}{
Results of running findBATCH to assess biological effect.
}
  \item{theme}{
ggplot theme
}
}
\details{
It searches for directions in the principal subspace that are associated with the batch variable.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A., 2017. A novel and robust statistical method to diagnose and correct batch effects in genomic data.
}
\author{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A.
}
