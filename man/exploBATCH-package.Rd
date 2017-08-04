\name{exploBATCH-package}
\alias{exploBATCH-package}
\alias{exploBATCH}
\docType{package}
\title{
\packageTitle{exploBATCH}
}
\description{
\packageDescription{exploBATCH}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{exploBATCH}
\packageIndices{exploBATCH}
Explore batch effect (exploBATCH) is an R package for perfoming formal statistical tests to discover batch effect and coorect 
for it. It searches for directions in the principal subspace that are associated with the batch variable of interest and correct for 
batch effect by subtracting this effect in the affected dimensions. Both the corrected and uncorrected dimensions are 
then used to predict the data.
}
\author{
\packageAuthor{exploBATCH}

Maintainer: \packageMaintainer{exploBATCH}
}
\references{
Nyamundanda, G., Poudel, P., Patil, Y. and Sadanandam, A., 2017. A novel and robust statistical method to diagnose and correct batch effects in genomic data.
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
}