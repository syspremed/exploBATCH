\name{expBATCH}
\alias{expBATCH}
\title{
Discovering and correcting batch effect.
}
\description{
A tool for discovering and correcting for batch effect using an approach in Nyamundanda et al (2017).
}
\usage{
expBATCH(D, batchCL, Conf = NA, mindim = 1, maxdim = 2, method = "ppcca", scale = "unit",SDselect=0)
}
\arguments{
  \item{D}{
Expression dataset in samples by features format.
}
  \item{batchCL}{
A vector identifying batches.
}
  \item{Conf}{
A vector identifying biological classes of interest.
}
  \item{mindim}{
Minimum number of principal components to be considered. Default is 1.
}
  \item{maxdim}{
Maximum number of principal components to be considered. Default is 2. maxdim should be greater or equal to mindim.
}
  \item{method}{
Method for batch correction. Either PPCCA or ComBat. PPCCA is the default.
}
\item{scale}{
Unit scale the data or none for no scaling. Unit scaling is the default 
}
\item{SDselect}{
Filter number of genes by setting SD value greater than zero using SDselect. SDselect=0 is the default.
}
}
\details{
It searches for directions in the principal subspace that are associated with batch effect and correct for batch effect by subtracting this 
effect in the affected dimensions. Both the corrected and uncorrected dimensions are then used to predict the data.
}
\value{
The output of this pipeline is saved in your working directory.
}
\references{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A., 2017. A novel and robust statistical method to diagnose and correct batch effects in genomic data.
}
\author{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A.
}
\examples{
# Expression data and batch effect variable is required for this function. The first example in the Nyamundanda et al (2017) manuscript 
# requires the gene expression data for breast cancer samples, named Breast (which is in the exploBATCHbreast package on github), and 
# variable identifying batches named batchBreast (in the exploBATCH package) 
require(devtools)
install_github("gnyamundanda/exploBATCHbreast") # install_github requires devtools package
require(exploBATCHbreast)
data(Breast)
data(batchBreast)
expBATCH(D=Breast,batchCL=batchBreast,Conf=NA,mindim=2,maxdim=9,method="ppcca",SDselect=2)
# The second example in the manuscript requires the gene expression data for colon cancer samples, named Colon (which is in the exploBATCHcolon 
# package on github), and variables identifying batches and tumor status (biological effect), named batchColon and bioCL respectively, in the 
# exploBATCH package. 
#require(devtools)
#install_github("gnyamundanda/exploBATCHcolon") # install_github requires devtools package
#require(exploBATCHcolon)
#data(Colon)
#data(batchColon)
#data(bioCL) 
#expBATCH(D=Colon,batchCL=batchColon,Conf=bioCL,mindim=2,maxdim=9,method="ppcca",SDselect=2)
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line